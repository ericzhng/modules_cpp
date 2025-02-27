
cmake_minimum_required(VERSION 3.30)
project(test_spdlog)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
get_filename_component(EXTERN_PATH "${EXTERN_PATH}" ABSOLUTE)

set(SPDLOG_DIR "${EXTERN_PATH}/spdlog" CACHE STRING "")


include(FetchContent)
FetchContent_Declare(
  cpptrace
  GIT_REPOSITORY https://github.com/jeremy-rifkin/cpptrace.git
  GIT_TAG        v0.8.2 # <HASH or TAG>
)
FetchContent_MakeAvailable(cpptrace)

set(CPPTRACE_DIR "${EXTERN_PATH}/cpptrace" CACHE STRING "")

set(INCLUDE_DIRS
	"${SPDLOG_DIR}/include"
	)

# include header directory
include_directories(${INCLUDE_DIRS})

# # Find Boost
# find_package(Boost REQUIRED COMPONENTS stacktrace_basic)

# set(CMAKE_PREFIX_PATH "${EXTERN_PATH}/spdlog/lib/cmake/spdlog" ${CMAKE_PREFIX_PATH})
# set(CMAKE_PREFIX_PATH "${EXTERN_PATH}/cpptrace/lib/cmake/cpptrace" ${CMAKE_PREFIX_PATH})
# find_package(spdlog CONFIG REQUIRED)
# find_package(cpptrace CONFIG REQUIRED)


# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

cmake_policy(SET CMP0167 NEW)

#define FMT_UNICODE 0
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")

# ensures that even in Release mode, the compiler will generate enough debugging information for cpptrace to work
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
# Some optimizations, like function inlining, may interfere with stack trace generation. You can disable specific optimizations in Release mode that affect debugging
# /Oy- disables frame pointer omission.
# /fp:precise ensures more precise floating-point calculations that can be important for stack trace generation.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Oy- /fp:precise")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O1")  # Moderate optimization

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GL")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

# add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)

# Create executable
add_executable(spdlogWrap main.cpp 
				spdlog_wrap.cpp
				spdlog_wrap.h)

# target_include_directories(spdlogWrap PRIVATE
    # "${CPPTRACE_DIR}/$<CONFIG>/include"
	# )
# target_link_directories(spdlogWrap PRIVATE 
    # "${CPPTRACE_DIR}/$<CONFIG>/lib"
	# )
	
target_link_libraries(spdlogWrap PRIVATE cpptrace::cpptrace dbghelp)

target_link_directories(spdlogWrap PRIVATE "${SPDLOG_DIR}/lib")
target_link_libraries(spdlogWrap PRIVATE $<IF:$<CONFIG:Debug>,spdlogd,spdlog>)

# target_link_libraries(spdlogWrap PRIVATE Boost::stacktrace_basic)
# # Enable backtrace-based stacktraces (optional)
# target_compile_definitions(spdlogWrap PRIVATE BOOST_STACKTRACE_ENABLE_BACKTRACE)

target_include_directories(spdlogWrap PRIVATE ${SPDLOG_INCLUDE_DIRS})
target_include_directories(spdlogWrap PRIVATE ${CPPTRACE_INCLUDE_DIRS})
