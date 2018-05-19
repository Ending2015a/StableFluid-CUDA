CXX = g++
LDFLAGS = -lm
CXXFLAGS = -O3 -std=c++11
TARGETS = main


.PHONY: all
all: $(TARGETS)

%: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -rf $(TARGETS)
