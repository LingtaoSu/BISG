NVCC := nvcc -std=c++11
CC := gcc -std=c99
CXX := g++ -std=c++11
LINK := g++
USEGPU = yes
DEBUG = no


# BLAS/LAPACK implementation
# Change these lines if you don't have an MKL but an other BLAS/LAPACK
LIBS=-llapack -lblas
#INCPATH=-I/opt/intel/mkl/include/
#LIBS=-L/opt/intel/mkl/lib/intel64/ -lmkl_rt

VPATH = librfn

GPUINC=-I/usr/local/cuda/include
GPULIB=-L/usr/local/cuda/lib64
GPUSO=-lcublas -lcurand -lcuda -lcudart -lcusolver -lgomp -lcusparse


ifeq ($(DEBUG), no)
	CFLAGS=-O3 -DNDEBUG -Wall -fPIC -march=native
	LDFLAGS=-O3 -flto -Wall -fPIC
else
	CFLAGS=-g -Wall -fPIC -march=native $(INCPATH)
	LDFLAGS= -g -Wall -fPIC $(LIBPATH) $(LIBS)
endif

ifeq ($(USEGPU),yes)
	INCPATH+=$(GPUINC)
	LIBS+=$(GPULIB) $(GPUSO)
else
	CFLAGS+=-DNOGPU
endif

CFLAGS+=$(INCPATH)
LDFLAGS+=$(LIBPATH) $(LIBS)
CXXFLAGS=$(CFLAGS)


# uncomment needed architectures as required
NVCCFLAGS=--use_fast_math $(addprefix -Xcompiler , $(CXXFLAGS)) \
           -gencode arch=compute_30,code=sm_35 \
           -gencode arch=compute_50,code=sm_50 \
           -gencode arch=compute_52,code=sm_52 
					 #-gencode arch=compute_61,code=sm_61


SOURCES=librfn.cpp cpu_operations.cpp nist_spblas.cc
OBJECTS=librfn.o cpu_operations.o nist_spblas.o

ifeq ($(USEGPU),yes)
	SOURCES+=gpu_operations.cu
	OBJECTS+=gpu_operations.o
endif

all: $(SOURCES) librfn.so

test: librfn.so tests/tests.o tests/test_runner.o
	g++ $(LDFLAGS) $^ -o $@ $(LIBS) -L./ -lrfn
	./test

testbin: librfn.so tests/testbin.o
	gcc tests/testbin.o -o testbin $(LIBPATH) $(LDFLAGS) -L./ -lrfn

test_rfn: librfn.so tests/test_rfn.o
	gcc tests/test_rfn.o -o test_rfn $(LIBPATH) $(LDFLAGS) -L./ -lrfn

librfn.so: $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS) -shared

gpu_operations.o: gpu_operations.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

clean:
	rm -rf *.o librfn.so tests/*.o
