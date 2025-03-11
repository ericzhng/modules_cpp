
cmake_minimum_required(VERSION 3.25)
project(MPI_laplace LANGUAGES Fortran)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_Fortran_INCLUDE_PATH})

# Create executable
add_executable(laplace_serial laplace_serial.f90)
target_link_libraries(laplace_serial   MPI::MPI_Fortran OpenMP::OpenMP_Fortran)

add_executable(laplace_mpi laplace_mpi.f90)
target_link_libraries(laplace_mpi      MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
