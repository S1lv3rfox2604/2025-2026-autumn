cmake_minimum_required(VERSION 3.10)

project(GaussSolver)

enable_language(CXX)

add_library(solver STATIC solver.cpp)