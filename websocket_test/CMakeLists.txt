
cmake_minimum_required(VERSION 3.30)
project(WebSocketServer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use vcpkg for package management
if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
endif()

# Find Boost WebSocket++
find_package(boost_asio REQUIRED CONFIG)
find_package(websocketpp CONFIG REQUIRED)

# Create executable
add_executable(WebSocketServer src/main.cpp src/WebSocketServer.cpp src/WebSocketServer.hpp)

# Link Boost and WebSocket++
target_include_directories(WebSocketServer PRIVATE ${Boost_INCLUDE_DIRS} ${websocketpp_INCLUDE_DIRS})
target_link_libraries(WebSocketServer PRIVATE websocketpp::websocketpp)
target_link_libraries(WebSocketServer PRIVATE Boost::asio)
