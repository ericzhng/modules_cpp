
cmake_minimum_required(VERSION 3.30)
project(python_c_link)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)

#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# # Use vcpkg for package management
# if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
#     message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
# endif()

set(CMAKE_TOOLCHAIN_FILE "C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake")

#============================================== #
# Find swig
# set(SWIG_DIR "C:/Portable/swigwin-4.0.2")
set(SWIG_EXECUTABLE "C:/Portable/swigwin-4.0.2/swig.exe")
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})


#============================================== #
# Find Python

set(Python_EXECUTABLE "D:/2-code/1-cpp-modules/extern/venv-python312/Scripts/python.exe" CACHE FILEPATH "Python executable")
set(Python_ROOT_DIR "C:/Python312")
set(Python_INCLUDE_DIR "${Python_ROOT_DIR}/include")
set(Python_LIBRARY "${Python_ROOT_DIR}/libs/python312.lib")
find_package(Python REQUIRED COMPONENTS Interpreter Development)

message(STATUS "Python Include: ${Python_INCLUDE_DIRS}")
message(STATUS "Python Lib: ${Python_LIBRARIES}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES mycpp.cpp) # Replace example.c with your C/C++ source files
set(SWIG_INTERFACE_FILE mycpp.i) # Replace example.i with your SWIG interface file

set_property(SOURCE mycpp.i PROPERTY CPLUSPLUS ON)
swig_add_library(mycpp 
				TYPE SHARED
				LANGUAGE python
				SOURCES ${SWIG_INTERFACE_FILE} ${SOURCE_FILES})
target_include_directories(mycpp PRIVATE "${Python_INCLUDE_DIRS}")
SWIG_LINK_LIBRARIES(mycpp ${Python_LIBRARIES})
# target_link_libraries(mycpp PRIVATE ${Python_LIBRARIES})
# set_target_properties(mycpp PROPERTIES PREFIX "" SUFFIX ".pyd") # For Windows, create a .pyd file

# for error message with python: cannot open file 'python312_d.lib'
# target_compile_definitions(mycpp PRIVATE SWIG_PYTHON_INTERPRETER_NO_DEBUG)

#============================================== #
# Another approach

set(MODULE_NAME "mycpp2")

# Custom command to run SWIG
add_custom_command( 
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/mycpp_wrap.cxx ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.py
    COMMAND ${SWIG_EXECUTABLE} -module ${MODULE_NAME} -python -c++ -o ${CMAKE_CURRENT_BINARY_DIR}/mycpp_wrap.cxx ${SWIG_INTERFACE_FILE}
    DEPENDS ${SWIG_INTERFACE_FILE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running SWIG to generate Python bindings"
    VERBATIM
)

# Create executable
add_library(mycpp_wrap MODULE ${CMAKE_CURRENT_BINARY_DIR}/mycpp_wrap.cxx mycpp.cpp)

set_target_properties(mycpp_wrap PROPERTIES PREFIX "" SUFFIX ".pyd" OUTPUT_NAME "_${MODULE_NAME}")  # Ensure correct .pyd extension on Windows

target_include_directories(mycpp_wrap PRIVATE "${Python_INCLUDE_DIRS}")
target_link_libraries(mycpp_wrap PRIVATE ${Python_LIBRARIES})

# for error message with python: cannot open file 'python312_d.lib'
# target_compile_definitions(mycpp_wrap PRIVATE SWIG_PYTHON_INTERPRETER_NO_DEBUG)
