Computing Homography Transformations
Christopher Nelson
11328750
christopher.w.nelso@email.wsu.edu

README.txt
LUdecomp.c
LUdecomp.h
homography.c
hmap.c
Makefile

This programming assignment uses the LUdecomposition engine from the last assignment to discover a transformation that maps the pixels from one image to another image. This transformation is called a homography and is performed in homography.c when N=4. To run, type make and then type:

./homography < input.in | ./hmap input.ppm xSize ySize > face.ppm 
