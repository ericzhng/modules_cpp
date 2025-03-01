
cmake_minimum_required(VERSION 3.30)
project(python_c_link)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)


#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()


#============================================== #
# Find packages
find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

pybind11_add_module(myadd myadd.cpp)
set_target_properties(myadd PROPERTIES PREFIX "" SUFFIX ".pyd")


#============================================== #
# Python method:

Python_add_library(myadd2 myadd2.cpp)
target_link_libraries(myadd2 PRIVATE pybind11::headers)
set_target_properties(myadd2 PROPERTIES
	INTERPROCEDURAL_OPTIMIZATION ON
	CXX_VISIBILITY_PRESET ON
	VISIBILITY_INLINES_HIDDEN ON
)
