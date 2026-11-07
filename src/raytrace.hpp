#pragma once
#include "geometry.hpp"
#include <cmath>

inline Vec3D refract(Vec3D incoming_ray, Vec3D normal, double eta, bool& internal_reflection){
	double c_1 = normal.dot(incoming_ray);
	double s_sqred = pow(eta, 2) * (1-pow(c_1, 2));
	if(s_sqred > 1.0){
		internal_reflection = true;
		return Vec3D{0,0,0};
	}
	double c_2 = std::sqrt(1-s_sqred);
	return (incoming_ray*eta) - normal*(eta * c_1 - c_2);
}
