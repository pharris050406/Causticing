#include "geometry.hpp"
#include "distributions.hpp"
#include "surface.hpp"
#include "raytrace.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

void check(const std::string& name, double got, double expected) {
    bool ok = std::fabs(got - expected) < 1e-12;
    std::cout << (ok ? "PASS " : "FAIL ") << name
              << ": got " << got << ", expected " << expected << '\n';
}

// k x k sites at the centers of a k x k grid, all weights 0.
// Every cell should be exactly its own small square, area 1/k^2.
void test_grid(int k) {
    std::vector<Vec2D> sites;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            sites.push_back({(i + 0.5) / k, (j + 0.5) / k});

    std::vector<Polygon> cells = power_diagram(sites, std::vector<double>(sites.size()));
    std::string name = "grid " + std::to_string(k) + "x" + std::to_string(k);

    check(name + " cell count", cells.size(), sites.size());

    // only print the cells that fail, so a 16x16 grid doesn't print 256 lines
    double worst = 0;
    for (std::size_t c = 0; c < cells.size(); c++) {
        double err = std::fabs(area(cells[c]) - 1.0 / (k * k));
        if (err > 1e-12) check(name + " cell " + std::to_string(c), area(cells[c]), 1.0 / (k * k));
        worst = std::fmax(worst, err);
    }
    check(name + " worst cell-area error", worst, 0.0);
}

// Random sites and weights: the cells must tile the square exactly,
// so their areas add up to 1 no matter how many cells come out empty.
void test_partition() {
    std::vector<Vec2D> sites = dist::uniform_random(500, 12345);
    std::vector<double> weights = dist::random_weights(sites.size(), 0.5, 999);
    std::vector<Polygon> cells = power_diagram(sites, weights);

    double total = 0;
    int empty = 0;
    for (const Polygon& cell : cells) {
        total += area(cell);
        if (cell.size() < 3) empty++;
    }
    std::cout << "partition: " << empty << " of " << cells.size() << " cells empty\n";
    check("partition total area", total, 1.0);
}

void test_surface() {
    const int k = 20;
    const double a = 1.0 / k, PI = 3.14159265358979323846;

    // 1. flat normals -> flat surface
    std::vector<Vec3D> flat(k * k, Vec3D{0, 0, 1});
    std::vector<double> h0 = integrate_heights(k, flat);
    double worst0 = 0;
    for (double h : h0) worst0 = std::fmax(worst0, std::fabs(h));
    check("surface flat", worst0, 0);

    // 2. constant slope p = 0.1 -> a plane rising 0.1*a per column
    std::vector<Vec3D> tilt(k * k, Vec3D{-0.1, 0, 1}.normalized());
    std::vector<double> h1 = integrate_heights(k, tilt);
	check("surface plane step", std::fabs((h1[1*k + 3] - h1[0*k + 3]) - 0.1*a) < 1e-9, 1);
    // 3. round trip: known surface -> exact normals -> integrate -> same surface?
    auto f = [&](double x, double y) { return 0.01 * std::sin(2*PI*x) * std::sin(2*PI*y); };
    std::vector<Vec3D> nrm; std::vector<double> exact;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            double x = (i + 0.5) * a, y = (j + 0.5) * a;
            double p = 0.01 * 2*PI * std::cos(2*PI*x) * std::sin(2*PI*y);
            double q = 0.01 * 2*PI * std::sin(2*PI*x) * std::cos(2*PI*y);
            nrm.push_back(Vec3D{-p, -q, 1}.normalized());
            exact.push_back(f(x, y));
        }
    std::vector<double> h2 = integrate_heights(k, nrm);
    double mean = 0; for (double e : exact) mean += e; mean /= exact.size();
    double worst2 = 0;
    for (std::size_t n = 0; n < h2.size(); n++) worst2 = std::fmax(worst2, std::fabs(h2[n] - (exact[n] - mean)));
    std::cout << "surface round trip: max error " << worst2 << " (surface amplitude 0.01)\n";
    check("surface round trip < 1e-4", worst2 < 1e-4, 1);
}
void test_raytrace() {
    const double eta = 1.49, D = 5;
    bool tir = false;
    const Vec3D up{0, 0, 1};

    // 1. flat surface: the ray passes straight through
    Vec3D t0 = refract(up, up, eta, tir);
    check("refract flat: no TIR", tir, 0);
    check("refract flat x", t0.x, 0); check("refract flat y", t0.y, 0); check("refract flat z", t0.z, 1);

    // 2. inverse of normal_to_target: the tilt it computes must send the ray to the target
    Vec3D n = normal_to_target({0.5, 0.5}, {0.7, 0.5}, D, eta);
    Vec3D t1 = refract(up, n, eta, tir);
    check("refract inverse: no TIR", tir, 0);
    check("refract inverse unit length", t1.length(), 1);
    check("refract inverse y", t1.y, 0);
    check("refract inverse aims at target", t1.x / t1.z, 0.2 / D);

    // 3. steep surface: the ray can't get out (total internal reflection)
    Vec3D steep = Vec3D{1, 0, 0.5}.normalized();
    refract(up, steep, eta, tir);
    check("refract steep: TIR", tir, 1);
}


int main() {
    // area
    check("unit square",      area({{0,0}, {1,0}, {1,1}, {0,1}}),  1.0);
    check("triangle",         area({{0,0}, {1,0}, {0,1}}),         0.5);
    check("empty",            area({}),                            0.0);
    check("clockwise square", area({{0,0}, {0,1}, {1,1}, {1,0}}), -1.0);
    check("shifted square",   area({{3,7}, {4,7}, {4,8}, {3,8}}),  1.0);

    // clip
    const Polygon sq = {{0,0}, {1,0}, {1,1}, {0,1}};
    check("clip x<=0.5",       area(clip(sq, {{ 1,0},  0.5})), 0.5);
    check("clip x>=0.5",       area(clip(sq, {{-1,0}, -0.5})), 0.5);
    check("clip thru corners", area(clip(sq, {{ 1,1},  1.0})), 0.5);
    check("clip corner tri",   area(clip(sq, {{ 1,1},  0.5})), 0.125);
    check("clip chop corner",  area(clip(sq, {{ 1,1},  1.5})), 0.875);
    check("clip all inside",   area(clip(sq, {{ 1,0},  2.0})), 1.0);
    check("clip all outside",  area(clip(sq, {{ 1,0}, -1.0})), 0.0);

    // power diagram
    std::vector<Vec2D> two = {{0.25,0.5}, {0.75,0.5}};
    auto c1 = power_diagram(two, {0.0, 0.0});
    check("2 sites equal w, cell 0", area(c1[0]), 0.5);
    auto c2 = power_diagram(two, {0.1, 0.0});
    check("2 sites w1-w2=0.1, cell 0", area(c2[0]), 0.6);
    check("2 sites w1-w2=0.1, cell 1", area(c2[1]), 0.4);
    auto c3 = power_diagram(two, {1.1, 1.0});
    check("2 sites shifted weights", area(c3[0]), 0.6);
    check("1 site", area(power_diagram({{0.3,0.8}}, {5.0})[0]), 1.0);

    test_grid(1);
    test_grid(4);
    test_grid(16);
    test_partition();
	check("dot", Vec3D{1,2,3}.dot(Vec3D{4,5,6}), 32);
	Vec3D cross1 = Vec3D{1,0,0}.cross(Vec3D{0,1,0});
	check("cross x", cross1.x, 0); check("cross y", cross1.y, 0); check("cross z", cross1.z, 1);

	Vec3D cross2 = Vec3D{0,1,0}.cross(Vec3D{1,0,0});
	check("cross reversed z", cross2.z, -1);
	check("length", Vec3D{3,4,12}.length(), 13);
	check("normalized z", Vec3D{0,0,5}.normalized().z, 1);
	// normal_to_target: screen at D = 5, acrylic eta = 1.49
	const double D = 5, eta = 1.49;

	Vec3D n0 = normal_to_target({0.5, 0.5}, {0.5, 0.5}, D, eta);      // no displacement: flat
	check("snell flat x", n0.x, 0); check("snell flat y", n0.y, 0); check("snell flat z", n0.z, 1);

	Vec3D np = normal_to_target({0.5, 0.5}, {0.7, 0.5}, D, eta);      // light sent toward +x
	check("snell unit length", np.length(), 1);
	check("snell +x: n_y is 0", np.y, 0);
	check("snell +x: n_x negative", np.x < 0, 1);                  // true counts as 1

	Vec3D nm = normal_to_target({0.5, 0.5}, {0.3, 0.5}, D, eta);      // mirror image: toward -x
	check("snell mirror x", nm.x, -np.x);
	check("snell mirror z", nm.z, np.z);

	// Snell's law itself: eta * sin(angle in) = sin(angle out), with sin = |n x d|
	Vec3D dI{0, 0, 1};
	Vec3D dT = Vec3D{0.2, 0, D}.normalized();
	check("snell's law", eta * np.cross(dI).length(), np.cross(dT).length());

	test_surface();
	test_raytrace();
}

