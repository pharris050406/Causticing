#pragma once

#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
struct DensityMap{
	int width, height;
	std::vector<double> densities;
};

inline DensityMap load_pgm(const std::string& path){
	std::ifstream file(path);
	if(!file){
		throw std::runtime_error("Unable to open file: cannot open " + path);
	}
	std::string magic;
	file >> magic;
	if(magic != "P2"){
		throw std::runtime_error("Unexpected file type; Expected P2 got \'" +magic + "\' instead.\nAre you using the right image (pgm)?");
	}
	int w, h;
	file >> w;
	file >> h;
	int maxval;
	file >> maxval;
	if(w<= 0 ||  h<= 0 || maxval <= 0){
		throw std::runtime_error("Unexpected file error: values/dimensions are negative");
	}

	DensityMap output{w, h, {}};
	int current_density;
	double sum = 0;
	// write pixel densities to DensityMap object
	for(int i = 0; i < (w*h); i++){
		file >> current_density;
		if(!file){
			throw std::runtime_error("File ended early");
		}
		output.densities.push_back(current_density);
		sum += current_density;
	}
	if(sum == 0){
		throw std::runtime_error("All black image, nothing to normalize.");
	}
	// normalize pixel values
	double mean = sum/(w*h);
	for(std::size_t i = 0; i < output.densities.size(); i++){
		output.densities.at(i) = output.densities.at(i)/mean;
	}
	return output;
}

inline void write_pgm(const std::string& path, const std::vector<double>& pixels, int width, int height){
	double brightest = 0;
	for(size_t i = 0; i < pixels.size(); i++){
		if(pixels[i] > brightest){
			brightest = pixels[i];
		}
	}

	std::ofstream file(path);
	if(!file){
		throw std::runtime_error("cannot open " + path + " for writing");
	}

	file << "P2\n" << width << " " << height << "\n255\n";
	for(int row = 0; row < height; row++){
		for(int col = 0; col < width; col++){
			double v = pixels[row * width + col];
			// scale so the brightest pixel is 255; an all-black image stays black
			int gray = (brightest > 0) ? static_cast<int>(std::lround(255.0 * v / brightest)) : 0;
			file << gray << (col + 1 < width ? ' ' : '\n');
		}
	}
}
