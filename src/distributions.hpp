#pragma once
// Site distributions and weight patterns to play with.
// Everything lives in the unit square [0,1]^2, same as power_diagram.
//
// Usage in main:
//   auto sites   = dist::sunflower(300);
//   auto weights = dist::wave_weights(sites, 0.01);
//   write_svg("sunflower.svg", power_diagram(sites, weights), sites);

#include "geometry.hpp"
#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <vector>

namespace dist {

constexpr double pi = std::numbers::pi;

// ---------------------------------------------------------------- sites

// n points uniformly at random. Change the seed for a different layout.
inline std::vector<Vec2D> uniform_random(int n, unsigned seed = 12345) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<Vec2D> s;
    for (int i = 0; i < n; ++i) s.push_back({u(rng), u(rng)});
    return s;
}

// k x k grid, each point nudged randomly inside its own grid square.
// jitter = 0 -> perfect squares, jitter = 1 -> points roam their whole square.
// Slide jitter from 0 to 1 to watch squares melt into a random diagram.
inline std::vector<Vec2D> jittered_grid(int k, double jitter, unsigned seed = 12345) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-0.5, 0.5);
    std::vector<Vec2D> s;
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            s.push_back({(i + 0.5 + jitter * u(rng)) / k,
                         (j + 0.5 + jitter * u(rng)) / k});
    return s;
}

// Hexagonal lattice: every other row shifted half a spacing -> honeycomb cells.
inline std::vector<Vec2D> hex_grid(int cols) {
    const double dx = 1.0 / cols;
    const double dy = dx * std::sqrt(3.0) / 2.0;
    std::vector<Vec2D> s;
    for (int row = 0; (row + 0.5) * dy < 1.0; ++row)
        for (int i = 0; i < cols; ++i)
            s.push_back({(i + 0.25 + 0.5 * (row % 2)) * dx, (row + 0.5) * dy});
    return s;
}

// Sunflower (Vogel) spiral: point i is turned by the golden angle from point i-1.
// Remarkably even spacing with no grid structure; the spiral arms show in the cells.
inline std::vector<Vec2D> sunflower(int n) {
    const double golden = pi * (3.0 - std::sqrt(5.0));
    std::vector<Vec2D> s;
    for (int i = 0; i < n; ++i) {
        double r = 0.48 * std::sqrt((i + 0.5) / n);
        double a = i * golden;
        s.push_back({0.5 + r * std::cos(a), 0.5 + r * std::sin(a)});
    }
    return s;
}

// Concentric rings around the center; ring r has points_inner * r points,
// so spacing along each ring stays roughly constant. Gives a dartboard.
inline std::vector<Vec2D> rings(int num_rings, int points_inner) {
    std::vector<Vec2D> s = {{0.5, 0.5}};
    for (int r = 1; r <= num_rings; ++r) {
        double radius = 0.45 * r / num_rings;
        int m = points_inner * r;
        for (int i = 0; i < m; ++i) {
            double a = 2.0 * pi * i / m + 0.3 * r;  // twist each ring a little
            s.push_back({0.5 + radius * std::cos(a), 0.5 + radius * std::sin(a)});
        }
    }
    return s;
}

// A few Gaussian blobs: tiny cells inside each cluster, huge cells between them.
// spread is the blob's standard deviation (try 0.03 - 0.1).
inline std::vector<Vec2D> clusters(int n, int num_clusters, double spread,
                                   unsigned seed = 12345) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.15, 0.85);
    std::normal_distribution<double> g(0.0, spread);
    std::vector<Vec2D> centers;
    for (int c = 0; c < num_clusters; ++c) centers.push_back({u(rng), u(rng)});
    std::vector<Vec2D> s;
    while (static_cast<int>(s.size()) < n) {
        const Vec2D& c = centers[s.size() % num_clusters];
        Vec2D p = {c.x + g(rng), c.y + g(rng)};
        if (p.x > 0.0 && p.x < 1.0 && p.y > 0.0 && p.y < 1.0) s.push_back(p);  // redraw if outside
    }
    return s;
}

// -------------------------------------------------------------- weights
// Weights have units of length^2, so "how big" depends on the site spacing.

// Random weights in [-strength, strength] / n  (1/n ~ squared spacing of n sites).
// strength 0.25: a few empty cells. 1: about a quarter empty. 2+: carnage.
inline std::vector<double> random_weights(std::size_t n, double strength,
                                          unsigned seed = 999) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> w(n);
    for (double& wi : w) wi = strength * u(rng) / n;
    return w;
}

// w = -2 (shift . s): every cell slides by `shift`, but the sites stay put.
// Exactly equivalent to the Voronoi diagram of the moved sites, so sites end up
// outside their own cells. Try shift = {0.1, 0} on a grid.
inline std::vector<double> slide_weights(const std::vector<Vec2D>& sites, Vec2D shift) {
    std::vector<double> w;
    for (const Vec2D& p : sites) w.push_back(-2.0 * shift.dot(p));
    return w;
}

// w = c |s - center|^2: the whole diagram is scaled by (1 - c) about the center.
// c = 0.3: interior cells shrink 30%, border cells grow to cover the rim.
// c < 0: the diagram expands and outer cells get pushed off the square (empty).
inline std::vector<double> bowl_weights(const std::vector<Vec2D>& sites, double c) {
    std::vector<double> w;
    for (const Vec2D& p : sites) {
        Vec2D d = p - Vec2D{0.5, 0.5};
        w.push_back(c * d.dot(d));
    }
    return w;
}

// Egg-crate pattern: alternating patches of big and small cells.
// amplitude is absolute: 0.003 is visible, 0.01 strong, 0.03 starts emptying cells.
// freq = number of waves across the square.
inline std::vector<double> wave_weights(const std::vector<Vec2D>& sites, double amplitude,
                                        double freq = 2.0) {
    std::vector<double> w;
    for (const Vec2D& p : sites)
        w.push_back(amplitude * std::sin(2 * pi * freq * p.x) * std::sin(2 * pi * freq * p.y));
    return w;
}

}  // namespace dist
