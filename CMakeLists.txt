cmake_minimum_required(VERSION 3.1)

project(mesh_shape_simplification)

set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    message(WARNING "OpenMP not found!")
endif()


set(SOURCES
        src/geom_utils.cpp
        src/geom_utils.h
        src/ply.h
        src/ply_utils.cpp
        src/ply_utils.h
        src/plyfile.c
)

add_executable(main src/main.cpp ${SOURCES})
