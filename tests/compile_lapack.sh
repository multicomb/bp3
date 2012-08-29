#!/bin/sh
g++ -O4 -g -funroll-loops -msse4 -o cholesky_inverse_lapack cholesky_inverse.cpp -llapack
#g++ -O4 -g -funroll-loops -mavx -o cholesky_inverse_lapack cholesky_inverse.cpp -llapack

