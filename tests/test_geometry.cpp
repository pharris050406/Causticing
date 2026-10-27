#include "geometry.hpp"
#include "distributions.hpp"
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
}
