
cmake_minimum_required(VERSION 3.30)
project(test_spdlog)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")

set(CMAKE_PREFIX_PATH "${EXTERN_PATH}/spdlog/lib/cmake/spdlog" ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH "${EXTERN_PATH}/cpptrace/lib/cmake/cpptrace" ${CMAKE_PREFIX_PATH})

set(INCLUDE_DIRS
	"${EXTERN_PATH}/spdlog/include"
	"${EXTERN_PATH}/cpptrace/include"
	)

# include header directory
include_directories(${INCLUDE_DIRS})

# Find Boost
find_package(Boost REQUIRED COMPONENTS stacktrace_basic)

find_package(spdlog CONFIG REQUIRED)
find_package(cpptrace CONFIG REQUIRED)

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

cmake_policy(SET CMP0167 NEW)

#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# Create executable
add_executable(spdlogWrap main.cpp 
				spdlog_wrap.cpp
				spdlog_wrap.h)

target_link_libraries(spdlogWrap PRIVATE Boost::stacktrace_basic)
target_link_libraries(spdlogWrap PRIVATE cpptrace::cpptrace)
target_link_libraries(spdlogWrap PRIVATE spdlog::spdlog)

# Enable backtrace-based stacktraces (optional)
target_compile_definitions(spdlogWrap PRIVATE BOOST_STACKTRACE_ENABLE_BACKTRACE)

# Link Boost and WebSocket++
target_include_directories(spdlogWrap PRIVATE ${SPDLOG_INCLUDE_DIRS})
target_include_directories(spdlogWrap PRIVATE ${CPPTRACE_INCLUDE_DIRS})
