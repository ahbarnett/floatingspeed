#!/bin/bash
# make and run some basic complex arithmetic and RAM access speed tests.
# Barnett 1/18/17, modified for Jeremy's tweaks.
# 2/3/17: Pataki suggested -Ofast.
# 3/28/17: realised that unrolling or march-native needed to get full speed for C++ complex vector std types on i7.

# clean
rm -f *.o complexmulttimingf complexmulttimingc complexmulttimingcpp

#export FLAGS="-O3"  # not enough, has bad (2x slower) hit in C++ complex types
#export FLAGS="-Ofast"

# either of following helps, esp for the final std C++ complex vector type...
#export FLAGS="-O3 -funroll-loops -march=native"

# actually slower than below:
#export FLAGS="-Ofast -funroll-loops -march=native"

# beats perf of Ofast, and allows correct Nan-checking, we hope:
#export FLAGS="-O3 -funroll-loops -march=native -fcx-limited-range"
export FLAGS="-O3 -march=native -fcx-limited-range"
# see discussion at: https://medium.com/@smcallis_71148/complex-arithmetic-is-complicated-873ec0c69fc5

#export FLAGS="-Ofast -funroll-loops -march=native"

# compile and run

GF=gfortran
#GF=gfortran-7

$GF complexmulttiming.f -o complexmulttimingf $FLAGS
echo Fortran:
./complexmulttimingf

# needed for complex.h if modern g++ version:
FLAGS+=" -std=c++03"

GXX=g++
#GXX=g++-7

g++-7 utils.cpp -c
$GXX complexmulttiming.cpp utils.o -o complexmulttimingc $FLAGS -D USE_C_TYPE_COMPLEX
echo C++-type real and C-type complex:
./complexmulttimingc

$GXX complexmulttiming.cpp utils.o -o complexmulttimingcpp $FLAGS -D USE_CPP_TYPE_COMPLEX
echo C++ type complex:
./complexmulttimingcpp
