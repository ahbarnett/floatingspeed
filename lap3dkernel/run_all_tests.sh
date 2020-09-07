#!/bin/bash
# compile and run all performance tests from linux. Barnett 3/4/20

echo C++SIMD:
# assumes g++ is compiler, but you'll want eg g++-9 for best single-thread:
(cd c++SIMD; g++ -fPIC -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./VCL_AgnerFog_version1+rsqrt -fopenmp lap3dkernel.cpp; ./a.out)

echo FORTRAN:
# assumes gfortran is your compiler; see above
(cd fortran; gfortran lap3dpottest.f90 -O3 -fopenmp -funroll-loops -march=native; ./a.out)


# uncomment following five lines to test intel c++ compiler and intel fortran compiler
#echo IC++SIMD:
#(cd c++SIMD; icpc -fPIC -O3 -march=native -funroll-loops -mkl -qopenmp -std=c++17 -DVCL -I./VCL_AgnerFog_version1+rsqrt -qopenmp lap3dkernel.cpp; ./a.out)
#
#echo IFORT:
#(cd fortran; ifort lap3dpottest.f90 -qopenmp -fPIC -O3 -march=native -funroll-loops -mkl; ./a.out)


echo JULIA with $(nproc) processors:
# note you'll need to import some Pkgs into your julia distro
JULIA_NUM_THREADS=$(nproc) julia --check-bounds=no -O3 julia/lap3dkernel.jl

echo PYTHON:
# assumes python is your v3 version w/ numba, eg via: source activate idp
python python/lap3dkernel.py

echo MATLAB:
# reminded how crappy its scripting behavior is (can't prevent startup text):
matlab -nojvm -nodisplay -log < matlab/lap3dkernel.m >/dev/null

echo OCTAVE:
octave matlab/lap3dkernel.m
