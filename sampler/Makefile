CXX = g++
LDFLAGS = -lm
CXXFLAGS = -O3 -std=c++11 -I./../lodepng -I./../cmdline
TARGETS = main
SOURCES = ../lodepng/lodepng.cpp

.PHONY: all
all: $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(SOURCES)

.PHONY: clean
clean:
	rm -rf $(TARGETS)
