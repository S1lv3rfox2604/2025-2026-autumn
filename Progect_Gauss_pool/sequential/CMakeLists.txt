cmake_minimum_required(VERSION 3.10)

project(SequentialSolver)

enable_language(CXX)

# БЕЗ оптимизаций - чистый код
add_library(sequential_solver STATIC solver.cpp)

# Устанавливаем выходной каталог
set_target_properties(sequential_solver PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)