cmake_minimum_required(VERSION 3.10)

project(Gauss)

enable_language(CXX)

add_subdirectory("init" ${CMAKE_CURRENT_BINARY_DIR}/initbin)
add_subdirectory("solver" ${CMAKE_CURRENT_BINARY_DIR}/solverbin)
add_subdirectory("print" ${CMAKE_CURRENT_BINARY_DIR}/printbin)

add_executable(exec main.cpp)

target_include_directories(exec PUBLIC "init")
target_link_libraries(exec PUBLIC init)

target_include_directories(exec PUBLIC "solver")
target_link_libraries(exec PUBLIC solver)

target_include_directories(exec PUBLIC "print")
target_link_libraries(exec PUBLIC print)

# Подключаем POSIX Threads
find_package(Threads REQUIRED)
target_link_libraries(exec PUBLIC Threads::Threads)
