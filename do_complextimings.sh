#!/bin/bash
# make and run some basic complex arithmetic and RAM access speed tests.
# Barnett 1/18/17, modified for Jeremy's tweaks.

# clean
rm -f *.o complexmulttimingf complexmulttiming

export FLAGS="-O3 -funroll-loops -march=native"
# compile and run
gfortran complexmulttiming.f -o complexmulttimingf $FLAGS
echo Fortran:
./complexmulttimingf
g++ utils.cpp -c
g++ complexmulttiming.cpp utils.o -o complexmulttiming $FLAGS -D USE_C_TYPE_COMPLEX
echo C and C-type:
./complexmulttiming
g++ complexmulttiming.cpp utils.o -o complexmulttiming $FLAGS -D USE_CPP_TYPE_COMPLEX
echo C++ type:
./complexmulttiming
