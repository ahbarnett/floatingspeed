#!/bin/bash
gfortran-9 lap3dpottest.f90 -O3 -fopenmp -funroll-loops -march=native
./a.out
