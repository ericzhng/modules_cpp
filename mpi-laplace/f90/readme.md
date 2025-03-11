
# use posix or linux to compile and test the code

cd "D:\2-code\1-cpp-modules\mpi-laplace\f90\build"

cd "D:\2-code\1-cpp-modules\mpi-pi-cal\build"

cmake -G "MinGW Makefiles" -DCMAKE_Fortran_COMPILER=mpifort ..

cmake --build .

## Not possibly to use it for fortran mpi in Windows using MSYS2
