#!/bin/bash
# make and run some basic complex arithmetic and RAM access speed tests.
# Barnett 1/18/17, modified for Jeremy's tweaks.
# 2/3/17: Pataki suggested -Ofast.

# clean
rm -f *.o complexmulttimingf complexmulttimingc complexmulttimingcpp

export FLAGS="-Ofast"
#export FLAGS="-Ofast -funroll-loops -march=native"
#export FLAGS="-O3"  # not enough

# compile and run

gfortran complexmulttiming.f -o complexmulttimingf $FLAGS
echo Fortran:
./complexmulttimingf

g++ utils.cpp -c
g++ complexmulttiming.cpp utils.o -o complexmulttimingc $FLAGS -D USE_C_TYPE_COMPLEX
echo C and C-type:
./complexmulttimingc

g++ complexmulttiming.cpp utils.o -o complexmulttimingcpp $FLAGS -D USE_CPP_TYPE_COMPLEX
echo C++ type:
./complexmulttimingcpp
