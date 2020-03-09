#!/bin/bash
# compile and run all performance tests from linux. Barnett 3/4/20

echo C++SIMD:
# assumes g++ is compiler, but you'll want eg g++-9 for best single-thread:
(cd c++SIMD; g++ -fPIC -O3 -march=native -funroll-loops -fopenmp -std=c++17 -DVCL -I./VCL_AgnerFog_version1 -fopenmp lap3dkernel.cpp; ./a.out)

echo FORTRAN:
# assumes gfortran is your compiler; see above
(cd fortran; gfortran lap3dpottest.f90 -O3 -fopenmp -funroll-loops -march=native; ./a.out)

echo JULIA:
# note you'll need to import some Pkgs into your julia distro
JULIA_NUM_THREADS=$(nproc) julia julia/lap3dkernel.jl

echo PYTHON:
# assumes python is your v3 version w/ numba, eg via: source activate idp
python python/lap3dkernel.py

echo MATLAB:
# reminded how crappy its scripting behavior is (can't prevent startup text):
matlab -nojvm -nodisplay -log < matlab/lap3dkernel.m >/dev/null

echo OCTAVE:
octave matlab/lap3dkernel.m
