
cmake_minimum_required(VERSION 3.30)
project(FEMSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# # Use vcpkg for package management
# if(CMAKE_TOOLCHAIN_FILE STREQUAL "")
    # message(FATAL_ERROR "Please define CMAKE_TOOLCHAIN_FILE to point to vcpkg.cmake")
# endif()

set(EXTERN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../extern")
set(EIGEN_INCLUDE_DIRS "${EXTERN_PATH}/eigen-3.4.0")  # Set this to your Eigen folder path

# set(CMAKE_PREFIX_PATH "${EXTERN_PATH}/eigen-3.4.0/cmake" ${CMAKE_PREFIX_PATH})
# find_package(Eigen3 REQUIRED)

# set(MPI_INCLUDE_DIR "C:/Program Files/Microsoft MPI/include")
# set(MPI_LIBRARIES "C:/Program Files/Microsoft MPI/lib/msmpi.lib")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MKL REQUIRED)

# set(MKL_INCLUDE_DIRS "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include")

# Add Eigen to your project (adjust path as necessary)
include_directories(${EIGEN_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")

# Create executable
add_executable(FEMSolver src/main.cpp src/FEMSolver.cpp src/FEMSolver.h)

# Link Boost and WebSocket++
target_include_directories(FEMSolver PRIVATE ${MKL_INCLUDE_DIRS})
target_link_libraries(FEMSolver ${MPI_LIBRARIES} MKL::MKL MPI::MPI_CXX OpenMP::OpenMP_CXX)

# --------------------------------------- #
# Create executable
add_executable(panua_test pardiso/test.cpp)

# Link Boost and WebSocket++
target_include_directories(panua_test PRIVATE "C:/tools/pardiso/include")
target_link_libraries(panua_test ${MPI_LIBRARIES} "C:/tools/pardiso/lib/libpardiso.lib" OpenMP::OpenMP_CXX)
add_custom_command(TARGET panua_test
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy "C:/tools/pardiso/lib/libpardiso.dll" ${CMAKE_BINARY_DIR}/.
)

# --------------------------------------- #
# Create executable
add_executable(pardiso_unsym_c pardiso/pardiso_unsym_c.cpp)
set_target_properties(pardiso_unsym_c PROPERTIES LINKER_LANGUAGE CXX)

# Link Boost and WebSocket++
target_include_directories(pardiso_unsym_c PRIVATE ${MKL_INCLUDE_DIRS})
target_link_libraries(pardiso_unsym_c ${MPI_LIBRARIES} MKL::MKL MPI::MPI_CXX OpenMP::OpenMP_CXX)

# --------------------------------------- #
# Create executable
add_executable(pardiso_sym_c pardiso/pardiso_sym_c.cpp)
set_target_properties(pardiso_sym_c PROPERTIES LINKER_LANGUAGE CXX)

# Link Boost and WebSocket++
target_include_directories(pardiso_sym_c PRIVATE ${MKL_INCLUDE_DIRS})
target_link_libraries(pardiso_sym_c ${MPI_LIBRARIES} MKL::MKL MPI::MPI_CXX OpenMP::OpenMP_CXX)
