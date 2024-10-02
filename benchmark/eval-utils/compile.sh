#! /bin/bash

g++ -std=c++11 -O3 set_intersection.cc -I/home/ylilo/projects/ust/PETS -o intersection -mavx2 -fopenmp