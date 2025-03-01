
cmake_minimum_required(VERSION 3.30)
project(python_c_link)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)

set(pybind11_DIR "${EXTERN_PATH}/pybind11")

#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# # Use vcpkg for package management
# if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
#     message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
# endif()


#============================================== #
# Find Python

set(Python_EXECUTABLE "D:/2-code/1-cpp-modules/extern/venv-python312/Scripts/python.exe" CACHE FILEPATH "Python executable")
set(Python_ROOT_DIR "C:/Python312")
set(Python_INCLUDE_DIR "${Python_ROOT_DIR}/include")
set(Python_LIBRARY "${Python_ROOT_DIR}/libs/python312.lib")
find_package(Python REQUIRED COMPONENTS Interpreter Development)

message(STATUS "Python Include: ${Python_INCLUDE_DIRS}")
message(STATUS "Python Lib: ${Python_LIBRARIES}")

# Create executable
add_library( myadd MODULE myadd.cpp )

set_target_properties(myadd PROPERTIES PREFIX "" SUFFIX ".pyd")

target_link_libraries(myadd PRIVATE "${Python_LIBRARIES}")
target_include_directories(myadd PRIVATE "${pybind11_DIR}/include" "${Python_INCLUDE_DIRS}")
