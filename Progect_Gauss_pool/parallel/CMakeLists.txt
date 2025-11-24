cmake_minimum_required(VERSION 3.10)

project(ParallelSolver)

enable_language(CXX)

# Флаги оптимизации
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -pthread")

add_library(parallel_solver STATIC solver.cpp)

# Устанавливаем выходной каталог
set_target_properties(parallel_solver PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)