#pragma once
#include "geometry.hpp"
#include <cstddef>
#include <raylib.h>
#include <vector>

inline Vector2 to_screen(const Vec2D& p){
	double x = p.x * GetScreenWidth();
	double y = (1-p.y) * GetScreenHeight();
	return Vector2{static_cast<float>(x), static_cast<float>(y)};
}
inline void draw_arrow(Vector2 from, Vector2 to, float head, Color color){
	float dx = to.x - from.x;
	float dy = to.y - from.y;
	float len = std::sqrt(dx*dx + dy*dy);
	DrawLineV(from, to, color);
	if(len < 1.0f){                  // too short to show a direction
		return;
	}
	dx /= len;  dy /= len;           // unit vector along the arrow
	float px = -dy, py = dx;         // unit vector perpendicular to it
	Vector2 left  = {to.x - head*dx + 0.5f*head*px, to.y - head*dy + 0.5f*head*py};
	Vector2 right = {to.x - head*dx - 0.5f*head*px, to.y - head*dy - 0.5f*head*py};
	DrawLineV(to, left, color);
	DrawLineV(to, right, color);
}
inline void draw_diagram(const std::vector<Polygon>& cells, const std::vector<Color>& colors){
	// pass 1: fill every cell
	for(std::size_t c = 0; c < cells.size(); c++){
		const Polygon& cell = cells[c];
		if(cell.size() < 3){
			continue;
		}
		// raylib only fills triangles that go counter-clockwise ON SCREEN.
		// our cells are CCW with y up; to_screen flips y, which makes them CW,
		// so hand the points over in reverse order.
		std::vector<Vector2> points;
		for(std::size_t i = 0; i < cell.size(); i++){
			points.push_back(to_screen(cell[i]));
		}
		DrawTriangleFan(points.data(), static_cast<int>(points.size()), colors[c]);
	}

	// pass 2: outlines on top
	for(const Polygon& cell : cells){
		if(cell.size() < 3){
			continue;
		}
		for(std::size_t i = 0; i < cell.size(); i++){
			DrawLineV(to_screen(cell[i]), to_screen(cell[(i+1) % cell.size()]), GRAY);
		}
	}
}
