#!/bin/bash
# make and run some basic complex arithmetic and RAM access speed tests.
# Barnett 1/18/17, modified for Jeremy's tweaks.

# clean
rm -f *.o complexmulttimingf complexmulttiming

# compile and run
gfortran complexmulttiming.f -o complexmulttimingf -O3 -funroll-loops
echo Fortran:
./complexmulttimingf
g++ utils.cpp -c
g++ complexmulttiming.cpp utils.o -o complexmulttiming -O3 -funroll-loops -D USE_C_TYPE_COMPLEX
echo C and C-type:
./complexmulttiming
g++ complexmulttiming.cpp utils.o -o complexmulttiming -O3 -funroll-loops -D USE_CPP_TYPE_COMPLEX
echo C++ type:
./complexmulttiming
