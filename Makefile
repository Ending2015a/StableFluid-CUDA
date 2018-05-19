NVCC = nvcc
LDFLAGS = -lm -lcublas -lcusparse
NVFLAGS = -O3 -std=c++11 -arch=sm_30
TARGETS = main
SOURCES = pcg_solver.cu

.PHONY: all
all: $(TARGETS)

%: %.cu
	$(NVCC) $(NVFLAGS) -Xcompiler="-O3" $(LDFLAGS) -o $@ $? $(SOURCES)

.PHONY: clean
clean:
	rm -rf $(TARGETS)
