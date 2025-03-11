
cmake_minimum_required(VERSION 3.30)
project(MPI_laplace LANGUAGES C)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

if (WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
elseif (UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


# no support of OpenACC on Windows
add_executable(laplace_mpi_acc c/laplace_mpi-grid_openacc_async_margin.c)
target_link_libraries(laplace_mpi_acc ${MPI_LIBRARIES})
