#pragma once
#include "geometry.hpp"
#include "image.hpp"
#include <cmath>
#include <algorithm>
#include <vector>

inline Vec3D refract(Vec3D incoming_ray, Vec3D normal, double eta, bool& internal_reflection){
	double c_1 = normal.dot(incoming_ray);
	double s_sqred = pow(eta, 2) * (1-pow(c_1, 2));
	if(s_sqred > 1.0){
		internal_reflection = true;
		return Vec3D{0,0,0};
	}
	internal_reflection = false;
	double c_2 = std::sqrt(1-s_sqred);
	return (incoming_ray*eta) - normal*(eta * c_1 - c_2);
}
// returns false if the point selected is outside the mesh
//
inline bool surface_at(int side_length, const std::vector<double>& heights, double x, double y, double& z, Vec3D& normal){
	double  a = 1.0/side_length;
	double u = (x - 0.5*a)/a;
	double v =(y-0.5*a)/a;
	if((u < 0 || u > side_length - 1) || (v < 0 || v > side_length - 1)){
		return false;
	}
	int i = std::min(static_cast<int>(std::floor(u)), side_length - 2);
	int j = std::min(static_cast<int>(std::floor(v)), side_length - 2);
	double fu = u - i;
	double fv = v - j;
	
	// find heights at each corner
	double h00 = heights[i * side_length + j];             // (i,   j)    bottom-left
	double h10 = heights[(i + 1) * side_length + j];       // (i+1, j)    bottom-right
	double h01 = heights[i * side_length + (j + 1)];       // (i,   j+1)  top-left
	double h11 = heights[(i + 1) * side_length + (j + 1)]; // (i+1, j+1)  top-right
	double p, q;
	if (fu + fv <= 1.0) {               // lower-left triangle: corners 00, 10, 01
		p = (h10 - h00) / a;
		q = (h01 - h00) / a;
		z = h00 + p * a * fu + q * a * fv;
	} else {                            // upper-right triangle: corners 11, 01, 10
		p = (h11 - h01) / a;
		q = (h11 - h10) / a;
		z = h11 - p * a * (1 - fu) - q * a * (1 - fv);
	}
	normal = Vec3D{-p, -q, 1}.normalized();
	return true;	
}


inline std::vector<double> trace_caustic(int side_length, const std::vector<double>& heights, double distance, double eta, int rays_per_side, int resolution){
	double a = 1.0/side_length; // site spacing
	double x_lo = 0.5 * a;
	double x_hi = 1 - 0.5 * a;
	double y_lo = 0.5 * a;
	double y_hi = 1 - 0.5 * a;

	std::vector<double> image = std::vector<double>(resolution*resolution, 0.0);
	int landed = 0;
	int off_screen = 0;
	int tir_count = 0;

	for(int r = 0; r < rays_per_side; r++){
		for(int c = 0; c < rays_per_side; c++){
			double x = x_lo + (r+0.5)*(x_hi-x_lo) / rays_per_side;
			double y = y_lo + (c+0.5)*(y_hi-y_lo) / rays_per_side;
			// these will get filled in on the next line
			double z; Vec3D n;
			if (!surface_at(side_length, heights, x, y, z, n)){
				continue;
			}
			bool tir;
			Vec3D t = refract(Vec3D{0,0,1}, n, eta, tir);
			if(tir){
				tir_count++;
				continue;
			}

			double s = (distance - z) / t.z; // how far along t until the ray reaches height `distance`
			double hx = x + s * t.x;
			double hy = y + s * t.y;
			if (hx >= 0 && hx < 1 && hy > 0 && hy <= 1) {
				int col = static_cast<int>(hx * resolution);
				int row = static_cast<int>((1 - hy) * resolution);   // flip: row 0 is the top
				image[row * resolution + col] += 1;
				landed++;
			} else {
				off_screen++;
			}
		}
	}
	std::cout << "rays landed " << landed << ", off screen " << off_screen
          << ", total internal reflection " << tir_count << '\n';
	return image;
}
