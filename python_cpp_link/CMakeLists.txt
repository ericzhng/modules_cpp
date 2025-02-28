
cmake_minimum_required(VERSION 3.30)
project(python_c_link)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)

# set(pybind11_DIR "${EXTERN_PATH}/pybind11")

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")


#============================================== #
# Create executable
add_library(ctypes SHARED ctypes/mycpp.cpp 
				)


#============================================== #
# Create executable
add_library(ffi SHARED ffi/mycpp.cpp 
				)


#============================================== #
# Find Boost
find_package(Boost REQUIRED COMPONENTS python)

# Create executable
add_library(boost_python SHARED boost-python/mycpp.cpp
				)
target_link_libraries(boost_python PRIVATE Boost::python)


#============================================== #
# Find Python
# set(Python_EXECUTABLE "C:/Software/anaconda3/envs/env_pipe/python.exe" CACHE FILEPATH "Python executable")
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
# Print paths for debugging
message(STATUS "Python Include: ${Python_INCLUDE_DIRS}")
message(STATUS "Python Lib: ${Python_LIBRARIES}")

# Create executable
add_library(swig_ex SHARED swig/mycpp_wrap.cxx
				)
include_directories("${Python_INCLUDE_DIRS}")
target_link_libraries(swig_ex PRIVATE ${Python_LIBRARIES})


#============================================== #
# Create executable
# Set paths to Python executable and libraries from your virtual environment
# set(Python_EXECUTABLE "C:/Software/anaconda3/envs/env_pipe/python.exe" CACHE FILEPATH "Python executable")
# set(Python_INCLUDE_DIR "C:/Software/anaconda3/envs/env_pipe/include" CACHE PATH "Python include directory")
# set(Python_LIBRARY "C:/Software/anaconda3/envs/env_pipe/libs/python312.lib" CACHE FILEPATH "Python library")

find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

# pybind11 method:
pybind11_add_module(pybind11_ex pybind11/mycpp.cpp)

# # Python method:
# Python_add_library(MyModule2 pybind11/mycpp.cpp)
# target_link_libraries(MyModule2 pybind11::headers)
# set_target_properties(MyModule2 PROPERTIES
	# INTERPROCEDURAL_OPTIMIZATION ON
	# CXX_VISIBILITY_PRESET ON
	# VISIBILITY_INLINES_HIDDEN ON
# )

