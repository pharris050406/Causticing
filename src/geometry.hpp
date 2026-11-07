#pragma once
#include <cmath>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <string>
#include <vector>

struct Vec2D {
  double x;
  double y;

  Vec2D operator+(const Vec2D &other) const {
    return {x + other.x, y + other.y};
  }
  Vec2D operator-(const Vec2D &other) const {
    return {x - other.x, y - other.y};
  }
  Vec2D operator*(double scalar) const { return {x * scalar, y * scalar}; }
  double dot(const Vec2D &vec) const {
	  return x * vec.x + y * vec.y;
  }

  std::string toString() const { return std::format("({}, {})", x, y); }
};

struct Vec3D{
	double x;
	double y;
	double z;
	
	Vec3D operator+(const Vec3D& other) const {
		return {x+other.x, y+other.y, z+other.z};
	}

	Vec3D operator-(const Vec3D &other) const {
		return {x - other.x, y - other.y, z - other.z};
	}
	Vec3D operator*(double scalar) const {
		return {x * scalar, y * scalar, z * scalar};
	}
	Vec3D operator/(double scalar) const {
		return {x / scalar, y / scalar, z / scalar};
	}
	Vec3D cross(const Vec3D& b) const {
		return Vec3D{y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x};
	};
	double dot(const Vec3D& b) const{
		return x*b.x + y*b.y + z*b.z;
	};
	double length() const{
		return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z,2));
	};
	Vec3D normalized(){
		return Vec3D{x,y,z}/Vec3D{x,y,z}.length();
	}
    std::string toString() const { return std::format("({}, {}, {})", x, y, z); }
};
// the refractive index is only of the material, since for air it's basically 1
inline Vec3D normal_to_target(Vec2D site, Vec2D target, double distance, double refIndex){
	Vec3D d_T = Vec3D{target.x - site.x, target.y - site.y, distance}.normalized();
	// vector representing the incoming direction of the ray through the glass.
	// aka up is through the glass
	Vec3D d_I = Vec3D{0,0,1}.normalized();
	return Vec3D (d_I*refIndex - d_T)/(d_I*refIndex - d_T).length();
}


struct HalfPlane {
  Vec2D normal;
  double offset;
};

inline Vec2D operator*(double scalar, Vec2D vector) {
  return (Vec2D{vector.x * scalar, vector.y * scalar});
}

using Polygon = std::vector<Vec2D>;
using Point = Vec2D;

inline double area(const Polygon &polygon) {
  double sum = 0;
  std::size_t size = polygon.size();
  if (size < 3) {
    return 0;
  }
  for (std::size_t i = 0; i < size; i++) {
    sum += (polygon.at(i).x * polygon.at((i + 1) % size).y -
            polygon.at((i + 1) % size).x * polygon.at(i).y);
  }

  return sum / 2;
}

inline Vec2D centroid(const Polygon& polygon){
	std::size_t num_vertices = polygon.size();
	if(num_vertices < 3){
		throw std::runtime_error("centroid of empty polygon");
	}
	double x = 0;
	double y = 0;
	for(size_t i = 0; i < num_vertices; i++){
		x+=((polygon[i].x + polygon[(i+1)%num_vertices].x)*(polygon[i].x * polygon[(i+1)%num_vertices].y - polygon[(i+1)%num_vertices].x * polygon[i].y));
		y+=((polygon[i].y + polygon[(i+1)%num_vertices].y)*(polygon[i].x * polygon[(i+1)%num_vertices].y - polygon[(i+1)%num_vertices].x * polygon[i].y));
	}
	x/=(area(polygon)*6);
	y/=(area(polygon)*6);
	return Vec2D{x,y};
}

inline Vec2D find_intersection(const Vec2D &start, const Vec2D &end,
                               const HalfPlane &halfPlane) {
  double t = (halfPlane.offset - halfPlane.normal.dot(start)) /
             (halfPlane.normal.dot(end - start));
  return start + t * (end - start);
}

inline Polygon clip(const Polygon &polygon, const HalfPlane &halfPlane) {
  Polygon output = {};
  std::size_t size = polygon.size();

  for (std::size_t i = 0; i < size; i++) {
    Vec2D start = polygon.at(i);
    Vec2D end = polygon.at((i + 1) % size);

    bool startInside = start.dot(halfPlane.normal) <= halfPlane.offset;
    bool endInside = end.dot(halfPlane.normal) <= halfPlane.offset;

    if (startInside && endInside) {
      output.push_back(end);
    } else if (startInside) {
      output.push_back(find_intersection(start, end, halfPlane));
    } else if (endInside) {
      output.push_back(find_intersection(start, end, halfPlane));
      output.push_back(end);
    }
  }
  return output;
}

inline std::vector<Polygon> power_diagram(const std::vector<Vec2D>& sites, const std::vector<double>& weights){
	std::vector<Polygon> output = {};
	std::size_t num_sites = sites.size();

	for(std::size_t i = 0; i < num_sites; i++){
		// begin each cell as a unit square. we're then going to just clip these down.
		Polygon cell = {{0,0},{1,0},{1,1},{0,1}};
		Vec2D site_i = sites.at(i);

		for(std::size_t j = 0; j < num_sites; j++){
			if(j==i){
				continue;
			}
			Vec2D site_j = sites.at(j);
			Vec2D normal = 2*(site_j-site_i);
			double offset = site_j.dot(site_j)-site_i.dot(site_i)-weights.at(j)+weights.at(i);
			HalfPlane halfPlane{normal, offset};
			cell = clip(cell, halfPlane);
			if(cell.size() < 3){
				break;
			}
		}
		output.push_back(cell);
	}

	return output;
}
