
cmake_minimum_required(VERSION 3.30)
project(ctypes-ffi)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# Create executable
add_library(mycpp SHARED mycpp.cpp)
