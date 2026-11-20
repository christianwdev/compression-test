# clang generates a noticeably faster decoder loop; g++ works too (make CXX=g++)
CXX := $(shell command -v clang++ >/dev/null 2>&1 && echo clang++ || echo g++)
CXXFLAGS ?= -O2 -g -std=c++17 -Wall -Wextra -Wno-misleading-indentation
SRC = src/main.cpp src/encoder.cpp src/decoder.cpp src/tables.cpp
build/gt: $(SRC) src/common.h
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -static -o $@ $(SRC)
clean:
	rm -rf build
