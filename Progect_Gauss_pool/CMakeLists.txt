cmake_minimum_required(VERSION 3.10)

project(Gauss)

enable_language(CXX)

# Используем встроенный поиск pthread вместо pkg-config
find_package(Threads REQUIRED)

# Создаем каталоги для выходных файлов
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Добавляем поддиректории в правильном порядке
add_subdirectory(sequential)
add_subdirectory(parallel)
add_subdirectory(init)
add_subdirectory(print)

add_executable(exec main.cpp)

# Подключаем sequential
target_include_directories(exec PUBLIC sequential)
target_link_libraries(exec sequential_solver)

# Подключаем parallel
target_include_directories(exec PUBLIC parallel)
target_link_libraries(exec parallel_solver)

# Подключаем остальные модули
target_include_directories(exec PUBLIC init)
target_link_libraries(exec init)

target_include_directories(exec PUBLIC print)
target_link_libraries(exec print)

# Связываем с библиотекой потоков
target_link_libraries(exec Threads::Threads)

# Добавляем флаги компиляции для pthread
target_compile_options(exec PRIVATE -pthread)