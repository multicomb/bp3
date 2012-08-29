#!/bin/sh
g++ -O4 -g -funroll-loops -msse4 -o cholesky_inverse_mkl cholesky_inverse.cpp -lmkl_core -lmkl_sequential -lmkl_intel_lp64 
# g++ -O4 -g -funroll-loops -mavx -o cholesky_inverse_mkl cholesky_inverse.cpp -lmkl_core -lmkl_sequential -lmkl_intel_lp64 

