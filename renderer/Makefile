CXX = g++
LDFLAGS = -lm -fopenmp
CXXFLAGS = -O3 -std=c++11
TARGETS = main
SOURCES = ../lodepng/lodepng.cpp

.PHONY: all
all: $(TARGETS)

%: %.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $< $(SOURCES)

.PHONY: clean
clean:
	rm -rf $(TARGETS)
