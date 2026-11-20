CXX ?= g++
CXXFLAGS ?= -O2 -g -std=c++17 -Wall -Wextra -Wno-misleading-indentation
SRC = src/main.cpp src/encoder.cpp src/decoder.cpp src/tables.cpp
build/gt: $(SRC) src/common.h
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)
clean:
	rm -rf build
