#pragma once
#include "geometry.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "image.hpp"

inline double mass_in_cell(const Polygon& cell, const DensityMap& img){
	if(cell.size() < 3){
		return 0.0;
	}
	double total = 0;
	double min_y = cell[0].y;
    double max_y = cell[0].y;
	double min_x = cell[0].x;
	double max_x = cell[0].x;
	for(Vec2D point:cell){
		if(point.x < min_x){
			min_x = point.x;
		}
		if(point.x > max_x){
			max_x = point.x;
		}
		if(point.y < min_y){
			min_y = point.y;
		}
		if(point.y > max_y){
			max_y = point.y;
		}
	}

	int px_lo = std::clamp(static_cast<int>(std::floor(min_x * img.width)),0,img.width - 1);
	int px_hi = std::clamp(static_cast<int>(std::floor(max_x * img.width)),0,img.width - 1);
	int py_hi = std::clamp(static_cast<int>(std::floor((1-min_y) * img.height)),0,img.height-1);
	int py_lo = std::clamp(static_cast<int>(std::floor((1-max_y) * img.height)),0,img.height-1);
	
	for(int py = py_lo; py <= py_hi; py++){
		for(int px = px_lo; px <= px_hi; px++){
			double	x0 = static_cast<double>(px)/img.width;
			double	x1 = static_cast<double>(px+1)/img.width;
			double	y0 = 1-static_cast<double>(py+1)/img.height;
			double	y1 = 1-static_cast<double>(py)/img.height;
			HalfPlane p1{Vec2D{-1,0}, -x0};
			HalfPlane p2{Vec2D{1,0}, x1};
			HalfPlane p3{Vec2D{0,-1}, -y0};
			HalfPlane p4{Vec2D{0,1}, y1};

			total+=area(clip(clip(clip(clip(cell, p1), p2), p3), p4)) * img.densities.at(py*img.width + px);

		}
	}
	return total;
}

inline Vec2D weighted_centroid(const Polygon& polygon, const DensityMap& img){
	size_t num_vertices = polygon.size();
	if(num_vertices < 3){
		throw std::runtime_error("polygon too small (< 3 vertices)");
	}
	double total_mass = 0;
	double centroid_x = 0;
	double centroid_y = 0;
	double min_y = polygon[0].y;
    double max_y = polygon[0].y;
	double min_x = polygon[0].x;
	double max_x = polygon[0].x;
	for(Vec2D point:polygon){
		if(point.x < min_x){
			min_x = point.x;
		}
		if(point.x > max_x){
			max_x = point.x;
		}
		if(point.y < min_y){
			min_y = point.y;
		}
		if(point.y > max_y){
			max_y = point.y;
		}
	}
	int px_lo = std::clamp(static_cast<int>(std::floor(min_x * img.width)),0,img.width - 1);
	int px_hi = std::clamp(static_cast<int>(std::floor(max_x * img.width)),0,img.width - 1);
	int py_hi = std::clamp(static_cast<int>(std::floor((1-min_y) * img.height)),0,img.height-1);
	int py_lo = std::clamp(static_cast<int>(std::floor((1-max_y) * img.height)),0,img.height-1);
	
	for(int py = py_lo; py <= py_hi; py++){
		for(int px = px_lo; px <= px_hi; px++){
			double	x0 = static_cast<double>(px)/img.width;
			double	x1 = static_cast<double>(px+1)/img.width;
			double	y0 = 1-static_cast<double>(py+1)/img.height;
			double	y1 = 1-static_cast<double>(py)/img.height;
			HalfPlane p1{Vec2D{-1,0}, -x0};
			HalfPlane p2{Vec2D{1,0}, x1};
			HalfPlane p3{Vec2D{0,-1}, -y0};
			HalfPlane p4{Vec2D{0,1}, y1};
			Polygon clipped = (clip(clip(clip(clip(polygon, p1), p2), p3), p4));
			if(clipped.size() < 3){
				continue;
			}
			double clipped_area = area(clipped);
			total_mass+=area(clipped) * img.densities.at(py*img.width + px);
			Vec2D c = centroid(clipped);
			centroid_x += img.densities[py*img.width + px] * clipped_area * c.x;
			centroid_y += img.densities[py*img.width + px] * clipped_area * c.y;
		}
	}
	if(total_mass <= 0){
		throw std::runtime_error("cell holds no light; no target");
	}

	return Vec2D{centroid_x / total_mass, centroid_y / total_mass};
}



// returns the Voronoi cell areas, one per site
// this is accurate for a uniform light distribution on a flat plane. 
inline std::vector<double> source_masses(const std::vector<Vec2D>& sites){
	std::vector<double> output;
	std::vector<Polygon> voronoi_diagram = power_diagram(sites, std::vector<double>(sites.size()));
	for(std::size_t i = 0; i < voronoi_diagram.size(); i++){
		output.push_back(area(voronoi_diagram.at(i)));
	}

	return output;
}
		
inline std::vector<double> mismatch(const std::vector<Vec2D>& sites, const std::vector<double>& weights, const DensityMap& image,const std::vector<double>& sources){
	std::vector<Polygon> diagram = power_diagram(sites, weights);
	// g calculates how much we're missing our target.
	// g[i] = light site i carries (sources[i]) - brightness its cell currently holds
	//   g > 0: cell holds too little -> raise w[i] to grow it
	//   g < 0: cell holds too much   -> lower w[i] to shrink it
	//   'sources' is fixed. it's the incoming light. 
	std::vector<double> g;
	for(std::size_t i = 0; i < diagram.size(); i++){
		g.push_back(sources[i] - mass_in_cell(diagram[i], image));
	}
	return g;
}
