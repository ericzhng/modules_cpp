
cmake_minimum_required(VERSION 3.30)
project(MPI_laplace LANGUAGES C)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

if (WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
elseif (UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


# Create executable
add_executable(laplace_serial laplace_serial.c utils_serial.h utils_timer.h)
target_link_libraries(laplace_serial MPI::MPI_C)

add_executable(laplace_mpi laplace_mpi.c utils_mpi.h utils_timer.h)
target_link_libraries(laplace_mpi MPI::MPI_C)
