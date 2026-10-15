#include <raylib.h>

int main(){
	InitWindow(800, 800, "caustics");
	while(!WindowShouldClose()){
		BeginDrawing();
			ClearBackground(RAYWHITE);
			DrawLineV(Vector2{0,0}, Vector2{100,100}, BLACK);
			DrawCircleV(Vector2{400,400}, 30, RED);
			EndDrawing();
	}
	CloseWindow();
	return 0;
}
