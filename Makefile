CXX      = c++
CXXFLAGS = -std=c++20 -Wall -Wextra -Isrc -I/usr/local/include
LDFLAGS  = -L/usr/local/lib -lraylib
DEBUG    = -g -fsanitize=address,undefined
HEADERS  = src/geometry.hpp src/image.hpp src/transport.hpp \
           src/distributions.hpp src/utils.hpp src/draw.hpp

run: build/relax
	./build/relax

test: build/test_geometry build/test_image
	./build/test_geometry
	./build/test_image

build/relax: src/main.cpp $(HEADERS)
	@mkdir -p build out
	$(CXX) $(CXXFLAGS) -O2 -o $@ src/main.cpp $(LDFLAGS)

build/test_geometry: tests/test_geometry.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(DEBUG) -o $@ tests/test_geometry.cpp

build/test_image: tests/test_image.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(DEBUG) -o $@ tests/test_image.cpp

clean:
	rm -rf build out

.PHONY: run test clean
