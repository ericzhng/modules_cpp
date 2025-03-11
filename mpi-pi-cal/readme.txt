
# how to compile, use cmd from msys2 ucrt64 environment

cmake -G "MinGW Makefiles" \
  -DMPI_C_COMPILER=/ucrt64/bin/mpicc \
  -DMPI_CXX_COMPILER=/ucrt64/bin/mpic++ \
  -DMPIEXEC="C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" \
  ..
  
cmake --build .
