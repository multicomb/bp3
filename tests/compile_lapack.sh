#!/bin/sh
g++ -O4 -g -funroll-loops -msse4 -Wall -o cholesky_inverse_lapack cholesky_inverse.cpp -llapack -fopenmp
#g++ -O4 -g -funroll-loops -mavx -Wall -o cholesky_inverse_lapack cholesky_inverse.cpp -llapack -fopenmp

