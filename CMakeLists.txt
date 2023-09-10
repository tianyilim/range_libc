cmake_minimum_required(VERSION 3.14)

PROJECT(range_lib
	VERSION 0.1
	DESCRIPTION "The ForzaETH fork of the RangeLibC ray-casting library"
	LANGUAGES CXX)

# Set CMake variables here
option(WITH_CUDA "Compile with CUDA ray cast variants" OFF)
# These variables are for CUDA
set(CHUNK_SIZE 262144) # this could be varied in the future
set(CHUNK_THREADS 256)

# Things that depend on top-level building go here
get_property(not_top DIRECTORY PROPERTY PARENT_DIRECTORY)

if(not_top)
	message(STATUS "RangeLibC building not in top-level build folder")

# TODO set equivalent output directories
else()
	message(STATUS "RangeLibC building in top-level build folder")
	set(CMAKE_BUILD_TYPE Release)

	# This should be set alongside the top-level colcon stuff, if possible
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
endif()

set_property(GLOBAL PROPERTY CXX_STANDARD 17)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED 17)

add_compile_options(-std=c++17 -Wall -Wextra -Wpedantic -march=native -fno-math-errno)

SET(BASEPATH "${CMAKE_SOURCE_DIR}")
add_definitions(-DBASEPATH=${BASEPATH})
add_subdirectory(${BASEPATH}/vendor/gflags)

# Use the basepath as a search location for includes
include_directories(${BASEPATH}/include ${BASEPATH}/vendor)

if(WITH_CUDA)
	add_definitions(-DUSE_CUDA=1)
	add_definitions(-DCHUNK_SIZE=${CHUNK_SIZE})
	add_definitions(-DNUM_THREADS=${CHUNK_THREADS})
	set(CUDA_CUDART_LIBRARY '/usr/local/cuda/lib64/libcudart.so')

	file(GLOB withCudaSrcs
		src/main.cpp
		"include/kernels.cu"
		"vendor/lodepng/lodepng.cpp")

	find_package(CUDA REQUIRED)

	cuda_add_executable(range_lib ${withCudaSrcs})
else()
	file(GLOB noCudaSrcs
		src/main.cpp
		"vendor/lodepng/lodepng.cpp")

	add_definitions(-DUSE_CUDA=0)
	add_executable(range_lib ${noCudaSrcs})
endif()

target_link_libraries(range_lib gflags pthread)