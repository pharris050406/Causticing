#pragma once
#include <cstddef>
#include <format>
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
  double dot(const Vec2D &vec) const { return x * vec.x + y * vec.y; }

  std::string toString() const { return std::format("({}, {})", x, y); }
};

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
