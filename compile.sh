#!/bin/bash

# setupATLAS
# lsetup "root 6.36.02-x86_64-el9-gcc14-opt"

COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
echo $COMPILER $FLAGS


$COMPILER $FLAGS -g -O3 -Wall -Wextra -Wpedantic -fopenmp ./draw1DHistos.cc  -Iinclude -o draw1DHistos
