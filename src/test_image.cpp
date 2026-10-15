
#include "image.hpp"
#include "transport.hpp"
#include "distributions.hpp"
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

void check(const std::string& name, double got, double expected) {
    bool ok = std::fabs(got - expected) < 1e-12;
    std::cout << (ok ? "PASS " : "FAIL ") << name
              << ": got " << got << ", expected " << expected << '\n';
}

// Passes if f() throws. Prints the exception's message so you can see it's the right one.
void expect_throw(const std::string& name, const std::function<void()>& f) {
    try {
        f();
        std::cout << "FAIL " << name << ": no exception thrown\n";
    } catch (const std::exception& e) {
        std::cout << "PASS " << name << ": threw \"" << e.what() << "\"\n";
    }
}

// Write a small PGM file to disk so each test controls the exact input.
void write_file(const std::string& path, const std::string& contents) {
    std::ofstream(path) << contents;
}

double mean_density(const Image& img) {
    double sum = 0;
    for (double d : img.densities) sum += d;
    return sum / (img.width * img.height);
}

// ---------------------------------------------------------------- mass_in_cell
// Images are built directly in code here, already normalized (mean density 1),
// so these tests don't depend on load_pgm.

void test_mass() {
    const Polygon square = {{0,0}, {1,0}, {1,1}, {0,1}};
    const Polygon left   = clip(square, {{ 1, 0},  0.5});   // x <= 0.5
    const Polygon right  = clip(square, {{-1, 0}, -0.5});   // x >= 0.5
    const Polygon top    = clip(square, {{ 0,-1}, -0.5});   // y >= 0.5
    const Polygon bottom = clip(square, {{ 0, 1},  0.5});   // y <= 0.5

    // uniform 4x4 image: mass should equal area for any cell
    Image uniform{4, 4, std::vector<double>(16, 1.0)};
    check("uniform: whole square",   mass_in_cell(square, uniform), 1.0);
    check("uniform: left half",      mass_in_cell(left, uniform), 0.5);
    check("uniform: corner triangle",
          mass_in_cell(clip(square, {{1, 1}, 0.5}), uniform), 0.125);
    check("uniform: off-grid sliver x<=0.3",                  // cuts through a pixel column
          mass_in_cell(clip(square, {{1, 0}, 0.3}), uniform), 0.3);

    // 2x1: left pixel white (density 2), right pixel black (density 0) -> tests x indexing
    Image lr{2, 1, {2.0, 0.0}};
    check("2x1: whole square", mass_in_cell(square, lr), 1.0);
    check("2x1: left half",    mass_in_cell(left, lr), 1.0);
    check("2x1: right half",   mass_in_cell(right, lr), 0.0);
    check("2x1: x<=0.25 (half of the white pixel)",
          mass_in_cell(clip(square, {{1, 0}, 0.25}), lr), 0.5);

    // 1x2: TOP row white, bottom row black -> tests the y flip (row 0 is the top)
    Image tb{1, 2, {2.0, 0.0}};
    check("1x2: top half",    mass_in_cell(top, tb), 1.0);
    check("1x2: bottom half", mass_in_cell(bottom, tb), 0.0);

    // partition: masses of all cells in a random power diagram sum to 1
    Image grad{8, 8, {}};                     // a smooth gradient, normalized below
    double s = 0;
    for (int i = 0; i < 64; i++) { grad.densities.push_back(i % 8 + 1); s += i % 8 + 1; }
    for (double& d : grad.densities) d /= s / 64;

    std::vector<Vec2D> sites = dist::uniform_random(200, 12345);
    std::vector<Polygon> cells = power_diagram(sites, dist::random_weights(sites.size(), 0.5, 999));
    double total = 0;
    for (const Polygon& c : cells) total += mass_in_cell(c, grad);
    check("partition: masses sum to 1 (gradient image)", total, 1.0);

    if (std::ifstream("images/target.pgm")) {
        Image t = load_pgm("images/target.pgm");
        double tt = 0;
        for (const Polygon& c : cells) tt += mass_in_cell(c, t);
        check("partition: masses sum to 1 (target.pgm)", tt, 1.0);
    }
}

int main() {
    try {
        // --- 2x1: one white pixel, one black pixel ---
        write_file("t_2x1.pgm", "P2\n2 1\n255\n255 0\n");
        Image a = load_pgm("t_2x1.pgm");
        check("2x1 width", a.width, 2);
        check("2x1 height", a.height, 1);
        check("2x1 pixel count", a.densities.size(), 2);
        check("2x1 white pixel density", a.densities.at(0), 2.0);
        check("2x1 black pixel density", a.densities.at(1), 0.0);

        // --- normalization on uneven values, spread over several lines ---
        write_file("t_3x2.pgm", "P2\n3 2\n255\n10 20 30\n40 50 60\n");
        Image b = load_pgm("t_3x2.pgm");
        check("3x2 pixel count", b.densities.size(), 6);
        check("3x2 mean density", mean_density(b), 1.0);
        check("3x2 first pixel", b.densities.at(0), 10.0 / 35.0);  // mean of raw values is 35

        // --- maxval doesn't matter after normalizing ---
        write_file("t_16bit.pgm", "P2\n2 1\n65535\n65535 0\n");
        check("16-bit white pixel density", load_pgm("t_16bit.pgm").densities.at(0), 2.0);

        // --- bad inputs must throw ---
        expect_throw("missing file",   [] { load_pgm("does_not_exist.pgm"); });
        write_file("t_p5.pgm", "P5\n2 1\n255\n");
        expect_throw("wrong magic",    [] { load_pgm("t_p5.pgm"); });
        write_file("t_zero.pgm", "P2\n0 1\n255\n");
        expect_throw("zero width",     [] { load_pgm("t_zero.pgm"); });
        write_file("t_trunc.pgm", "P2\n3 2\n255\n10 20 30\n40\n");
        expect_throw("truncated file", [] { load_pgm("t_trunc.pgm"); });
        write_file("t_black.pgm", "P2\n2 2\n255\n0 0\n0 0\n");
        expect_throw("all-black image", [] { load_pgm("t_black.pgm"); });

        // --- your real image, if it's in this directory ---
        if (std::ifstream("target.pgm")) {
            Image t = load_pgm("target.pgm");
            std::cout << "target.pgm: " << t.width << "x" << t.height << '\n';
            check("target.pgm pixel count", t.densities.size(), t.width * t.height);
            check("target.pgm mean density", mean_density(t), 1.0);
        } else {
            std::cout << "(no target.pgm here, skipping)\n";
        }

        test_mass();
    } catch (const std::exception& e) {
        std::cerr << "unexpected exception: " << e.what() << '\n';
        return 1;
    }
}
