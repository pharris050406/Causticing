#pragma once
#include "geometry.hpp"
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
inline void write_svg(const std::string& path, const std::vector<Polygon>& cells, const std::vector<Vec2D>& sites){
	std::ofstream outFile(path);
	if(!outFile){
		std::perror(path.c_str());
		return;
	}
	outFile << "<svg xmlns='http://www.w3.org/2000/svg' width='800' height='800'>\n";
	outFile << "<rect width='800' height='800' fill='white'/>\n";
	for(const Polygon& cell:cells){
		if(cell.size() < 3){
			continue;
		}
		outFile << "<polygon points='"; 
		for(const Vec2D& vertex:cell){
			outFile << (std::to_string(vertex.x * 800) + "," + std::to_string((1 - vertex.y) * 800) + " ");
		}
		outFile << "' fill='none' stroke='black'/>\n";
	}

	for(const Vec2D& site:sites){
		
		outFile << "<circle cx='";
		outFile << (std::to_string(site.x * 800) + "' cy='" + std::to_string((1 - site.y) * 800) + "' r='3'");
		outFile << " fill='red' stroke='none'/>\n";
	}
	outFile << "</svg>";
	outFile.close();
}

inline double sum_abs(const std::vector<double>& vec){
	double sum = 0;
	for(double item:vec){
		sum+=std::fabs(item);
	}
	return sum;
}

inline double max_abs(const std::vector<double>& vec){
	double highest = 0;
	for(double item:vec){
		if(std::fabs(item) > highest){
			highest = std::fabs(item);
		}
	}
	return highest;
}
