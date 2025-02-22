
cmake_minimum_required(VERSION 3.30)
project(NamedPipeEx)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# # Find Boost
# find_package(boost_asio REQUIRED CONFIG)

# Create executable
add_executable(NamedPipeEx main_named_pipe.cpp PipeWriter.cpp PipeWriter.h)
target_compile_options(NamedPipeEx PRIVATE -D_UNICODE) 

# # Link Boost and WebSocket++
# target_include_directories(NamedPipeEx PRIVATE ${Boost_INCLUDE_DIRS} ${NamedPipeEx_INCLUDE_DIRS})
# target_link_libraries(NamedPipeEx PRIVATE NamedPipeEx::NamedPipeEx)
# target_link_libraries(NamedPipeEx PRIVATE Boost::asio)
