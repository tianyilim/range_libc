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

add_compile_options(-std=c++17 -Wall -Wextra -Wpedantic -march=native -fno-math-errno -O2 -fPIC)

add_subdirectory(vendor/gflags)
add_subdirectory(vendor/pybind11)

# You can change this path to point where eigen3 lives on your system
include_directories(include vendor /usr/include/eigen3)

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

# This is a test driver
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