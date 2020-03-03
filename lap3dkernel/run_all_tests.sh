#!/bin/bash
# compile and run all performance tests from linux. Barnett 3/3/20

echo FORTRAN:
# assumes gfortran is your compiler; use module load GCC/7.4.0 etc
(cd fortran; gfortran lap3dpottest.f90 -O3 -fopenmp -funroll-loops -march=native; ./a.out)

echo JULIA:
JULIA_NUM_THREADS=$(nproc) julia julia/lap3dkernel.jl

echo PYTHON:
# assumes python is your v3 version w/ numba, eg via: source activate idp
python python/lap3dkernel.py

echo MATLAB:
# reminded how crappy its scripting behavior is (can't prevent startup text):
matlab -nojvm -nodisplay -log < matlab/lap3dkernel.m >/dev/null

echo OCTAVE:
octave matlab/lap3dkernel.m
