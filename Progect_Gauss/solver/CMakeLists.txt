cmake_minimum_required(VERSION 3.10)

project(GaussSolver)

enable_language(CXX)

add_library(solver STATIC solver.cpp)

# Подключаем POSIX Threads
find_package(Threads REQUIRED)
target_link_libraries(solver PUBLIC Threads::Threads)
