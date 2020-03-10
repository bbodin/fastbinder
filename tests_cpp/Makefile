C_FILES := $(shell find ../fastbinder/src -name \*.\*pp -not -name RcppExports.cpp)

R_FLAGS= -m64 -std=gnu++11  -DNDEBUG   -fopenmp   -O3 -g -pipe -Wall \
          -fexceptions -fstack-protector-strong -grecord-gcc-switches \
          -march=native -mtune=native -fasynchronous-unwind-tables  -lboost_program_options -lboost_system

I_FLAGS= -I/usr/local/cuda-10.0/include/  -I../fastbinder/src/ -I /usr/include/R/
L_FLAGS= -L/usr/local/cuda-10.0/lib64/ -lOpenCL

TEST_FILES=$(shell ls *.cpp)

TARGETS := test_mixtures test_binders test_mixtures_debug
all : ${TARGETS}

verbose.o :  ../fastbinder/src/verbose.cpp ${C_FILES} 
	g++ ${R_FLAGS}  ../fastbinder/src/verbose.cpp -DNO_RCPP  ${I_FLAGS}  -larmadillo  -c -o $@ 

binders.o :  ../fastbinder/src/binders.cpp ${C_FILES} 
	g++ ${R_FLAGS}   ../fastbinder/src/binders.cpp -DNO_RCPP  ${I_FLAGS}  -larmadillo  -c -o $@ 


test_binders : test_binders.cpp binders.o ./binders/opencl_binders.cpp  verbose.o
	g++ ${R_FLAGS} ${I_FLAGS}  ${L_FLAGS} $^ -DNO_RCPP  -larmadillo  -o $@ 


clean :
	rm ${TARGETS} *.o  -f 

