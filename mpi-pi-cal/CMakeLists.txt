
cmake_minimum_required(VERSION 3.30)



project(calc_pi_c LANGUAGES C)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# Create executable
add_executable(calc_c calc_pi_mpi.c)
target_link_libraries(calc_c MPI::MPI_C OpenMP::OpenMP_C)




project(calc_pi_f LANGUAGES Fortran)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# Create executable
add_executable(calc_f calc_pi_mpi.f90)
target_link_libraries(calc_f MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
