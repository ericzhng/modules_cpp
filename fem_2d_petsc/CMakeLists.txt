
cmake_minimum_required(VERSION 3.25)
project(FEMSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable optimization flags (optional)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")


set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
set(EIGEN_INCLUDE_DIRS "${EXTERN_PATH}/eigen-3.4.0")  # Set this to your Eigen folder path

# Add Eigen to your project (adjust path as necessary)
include_directories(${EIGEN_INCLUDE_DIRS})


# Find MPI (required by PETSc and your code)
find_package(MPI REQUIRED)
if(MPI_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_definitions(-DOMPI_SKIP_MPICXX) # Avoid deprecated MPI C++ bindings
else()
    message(FATAL_ERROR "MPI not found. Please install an MPI implementation.")
endif()


find_package(OpenMP REQUIRED)


# Find PETSc
set(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "Path to PETSc installation")
set(PETSC_ARCH $ENV{PETSC_ARCH} CACHE STRING "PETSc architecture (e.g., arch-linux-c-debug)")
if(NOT PETSC_DIR)
    message(FATAL_ERROR "PETSC_DIR not set. Please set the environment variable PETSC_DIR.")
endif()

# Include PETSc's configuration
# include(${PETSC_DIR}/lib/petsc/conf/petscvariables)
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)
if(PETSC_FOUND)
    include_directories(${PETSC_INCLUDE_DIRS})
    link_directories(${PETSC_LIBRARY_DIRS})
    message(STATUS "PETSc found: ${PETSC_LIBRARIES}")
else()
    message(FATAL_ERROR "PETSc not found. Check PETSC_DIR and PETSC_ARCH.")
endif()


# Create executable
add_executable(FEMSolver src/main.cpp src/FEMSolver.cpp src/FEMSolver.h)

# Link Boost and WebSocket++
target_link_libraries(FEMSolver ${PETSC_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)

