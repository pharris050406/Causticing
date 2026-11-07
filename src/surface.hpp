#pragma once
#include <cmath>
#include "geometry.hpp"
#include <vector>

inline std::vector<double> integrate_heights(int side_length, const std::vector<Vec3D>& normals){
	std::vector<double> heights(pow(side_length, 2), 0.0);
	std::vector<double> rise_x, rise_y;
	double spacing = 1.0/side_length;
	for(size_t i = 0; i < normals.size(); i++){
		rise_x.push_back((-1) * normals[i].x / normals[i].z);
		rise_y.push_back((-1) * normals[i].y / normals[i].z);
	}
	

	
	for(int sweep = 0; sweep < 10000; sweep++){
		double max_change = 0;
		for(int i = 0; i < side_length; i++){
			for(int j = 0; j < side_length; j++){
				int n = side_length;
				int current = i * n + j;
				int left = (i-1) * n + j;
				int right  = (i+1) * n + j;
				int below = i * n + (j-1);
				int above = i * n + (j+1);
				double sum = 0;
				int count = 0;

				if (i > 0) {        // left neighbor: rise from it to current, going right
					sum += heights[left] + spacing * (rise_x[left] + rise_x[current]) / 2;
					count++;
				}
				if (i < n - 1) {    // right neighbor: current is BELOW it by the rise
					sum += heights[right] - spacing * (rise_x[current] + rise_x[right]) / 2;
					count++;
				}
				if (j > 0) {        // below neighbor
					sum += heights[below] + spacing * (rise_y[below] + rise_y[current]) / 2;
					count++;
				}
				if (j < n - 1) {    // above neighbor
					sum += heights[above] - spacing * (rise_y[current] + rise_y[above]) / 2;
					count++;
				}

				double new_height = sum / count;
				double change = std::fabs(new_height - heights[current]);
				if(change > max_change){
					max_change = change;
				}

				heights[current] = new_height;
			}
		}
		if(max_change < 1e-12){
			break;
		}
	}
	double mean = 0.0;
	for(int i = 0; i < pow(side_length, 2); i++){
		mean += heights[i];
	}
	mean /= pow(side_length,2);

	for(double& h : heights){
		h-=mean;
	}

	return heights;
}

