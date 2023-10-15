cmake_minimum_required(VERSION 3.14)

PROJECT(range_lib_ray_casting
	VERSION 0.1
	DESCRIPTION "The ForzaETH fork of the RangeLibC ray-casting library"
	LANGUAGES CXX)

# GoogleTest min C++ version is 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set CMake variables here
option(WITH_CUDA "Compile with CUDA ray cast variants" OFF)

# These variables are for CUDA
set(CHUNK_SIZE 262144) # this could be varied in the future
set(CHUNK_THREADS 256)

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
	googletest
	URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
FetchContent_MakeAvailable(googletest)

set(CMAKE_BUILD_TYPE Release)

# Things that depend on top-level building go here
get_property(not_top DIRECTORY PROPERTY PARENT_DIRECTORY)

# In the end not needed, Rangelibc was able to build @ top level
if(not_top)
	message(STATUS "RangeLibC building not in top-level build folder")
else()
	message(STATUS "RangeLibC building in top-level build folder")

	# This should be set alongside the top-level colcon stuff, if possible
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
endif()

add_compile_options(-std=c++17 -Wall -Wextra -Wpedantic -march=native -fno-math-errno -O2 -fPIC)

add_subdirectory(vendor/gflags)
add_subdirectory(vendor/pybind11)

include_directories(include vendor)

add_library(lodepng_lib SHARED
	vendor/lodepng/lodepng.cpp
)

add_library(rangeUtils SHARED
	src/RangeUtils.cpp
)

add_library(geometryUtils SHARED
	src/geometry.cpp
)

pybind11_add_module(rangelib MODULE pywrapper/rangelib_pybindings.cpp)
target_link_libraries(rangelib PUBLIC lodepng_lib rangeUtils geometryUtils)

# if(WITH_CUDA)
# add_definitions(-DUSE_CUDA=1)
# add_definitions(-DCHUNK_SIZE=${CHUNK_SIZE})
# add_definitions(-DNUM_THREADS=${CHUNK_THREADS})
# set(CUDA_CUDART_LIBRARY '/usr/local/cuda/lib64/libcudart.so')

# file(GLOB withCudaSrcs
# src/main.cpp
# "include/kernels.cu")

# find_package(CUDA REQUIRED)

# cuda_add_executable(range_lib ${withCudaSrcs})
# else()
# file(GLOB noCudaSrcs
# src/main.cpp)

# add_definitions(-DUSE_CUDA=0)
# add_executable(range_lib ${noCudaSrcs})
# endif()

# target_link_libraries(range_lib gflags pthread lodepng_lib)
add_executable(rg_test
	src/testDriver.cpp
)
target_link_libraries(rg_test lodepng_lib rangeUtils geometryUtils)

# Testing only
enable_testing()

add_executable(
	test_rangelib

	test/test_omap.cpp
	test/test_bl.cpp
	test/test_glt.cpp
	test/test_ray_marching.cpp
	test/test_geom_utils.cpp
)
target_link_libraries(
	test_rangelib
	lodepng_lib
	rangeUtils
	geometryUtils
	GTest::gtest_main
)
include(GoogleTest)