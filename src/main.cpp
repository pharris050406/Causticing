#include "geometry.hpp"
#include "raylib.h"
#include "utils.hpp"
#include "distributions.hpp"
#include "transport.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "image.hpp"
#include "draw.hpp"
int main(){
	try{
	double step = 0.02;
	int num_sites = 600;
	bool stuck = false;

	InitWindow(800, 800, "caustics");
	SetTargetFPS(60);

	DensityMap main_image = load_pgm("images/cat.pgm");

	std::vector<Vec2D> sites = dist::jittered_grid(static_cast<int>(sqrt(num_sites)), 0);
	std::vector<double> weights(sites.size());
	std::vector<double> sources = source_masses(sites);
	std::vector<double> g = mismatch(sites, weights, main_image, sources);

	double d_max = *std::max_element(main_image.densities.begin(), main_image.densities.end());

	bool done = false;
	std::vector<Vec2D> targets;

	while(!WindowShouldClose()){
		if(!done){
			if(max_abs(g) * sites.size() < 0.01){
				std::cout << "image has converged\n";
				done = true;
			} else {
				while(true){
					std::vector<double> weights2 = weights;
					for(std::size_t j = 0; j < weights2.size(); j++){
						weights2[j] += step * g[j];
					}
					std::vector<double> g2 = mismatch(sites, weights2, main_image, sources);

					if(sum_abs(g2) < sum_abs(g)){
						weights = weights2;
						g = g2;
						step = step * 1.2;
						break;
					}else{
						step = step * 0.5;
					}

					if(step <= 1E-10){
						stuck = true;
						break;
					}
				}

				std::cout << "step: " + std::to_string(step) + "\n";
				std::cout << "max_abs(g) * N: " + std::to_string(max_abs(g) * sites.size()) + "\n";
				if(stuck){
					std::cout << "stuck!\n";
					done = true;
				}
			}
		}

		//drawing the frames (runs every frame, even after done)
		std::vector<Polygon> drawing_cells = power_diagram(sites, weights);

		// color each cell by the brightness it shows: light it carries / area it covers
		std::vector<Color> colors;
		for (std::size_t i = 0; i < drawing_cells.size(); i++) {
			double a = area(drawing_cells[i]);
			if (a <= 0) {                  // empty cell: never drawn, so its color doesn't matter
				colors.push_back(BLACK);
				continue;
			}
			double d = sources[i] / a;
			double v = 255.0 * std::min(1.0, d / d_max);   // clamp: unconverged cells can exceed d_max
			unsigned char c = static_cast<unsigned char>(v);
			colors.push_back(Color{c, c, c, 255});
		}

		if(done && targets.empty()){
			for(size_t i = 0; i < drawing_cells.size(); i++){
				if(drawing_cells[i].size() < 3){
					targets.push_back(sites[i]);
				}else{
					targets.push_back(weighted_centroid(drawing_cells[i], main_image));
				}
			}
		}

		BeginDrawing();
			ClearBackground(RAYWHITE);
			draw_diagram(drawing_cells, colors);
			for(size_t i = 0; i < targets.size(); i++){
				draw_arrow(to_screen(sites[i]), to_screen(targets[i]), 5, RED);
			}
				
		EndDrawing();
	}

	CloseWindow();

	std::vector<Polygon> cells = power_diagram(sites, weights);
	write_svg("out/new.svg", cells, sites);

	return 0;
	}catch(std::exception& e){
		std::cerr << e.what() << "\n";
		return 1;
	}
}
