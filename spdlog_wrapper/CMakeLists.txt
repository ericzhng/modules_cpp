
cmake_minimum_required(VERSION 3.30)
project(test_spdlog)

set(SPDLOG_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/spdlog-1.15.1/include")

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

# Create executable
add_executable(spdlogWrap spdlog_wrapper.cpp Logger.h)

# Link Boost and WebSocket++
target_include_directories(spdlogWrap PRIVATE ${SPDLOG_INCLUDE_DIRS})
