# Compiler and compilation flags
cc=gcc
cxx=g++
cflags=-Wall -ansi -pedantic -fopenmp -O3

# Compiler for ImageMagick-dependent executables. These must be defined separately
# due to linking issues on the Mac.
im_cxx=$(cxx)
#im_cflags=$(cflags)
#im_lflags=

#im_cflags=-DHAS_IMAGEMAGICK `Magick++-config --cppflags --cxxflags`
#im_lflags=`Magick++-config --ldflags --libs`

# MPI compiler
mpicxx=mpicxx -Wno-long-long

# LAPACK flags for dense linear algebra
lp_lflags=-llapack -lblas
