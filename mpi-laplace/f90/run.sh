
cmake -G "MinGW Makefiles" -B build \
  -DMPI_C_COMPILER=/ucrt64/bin/gcc \
  -DMPI_CXX_COMPILER=/ucrt64/bin/g++ \
  -DMPIEXEC="C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" \
  -DMPI_C_INCLUDE_PATH=/ucrt64/include \
  -DMPI_C_LIBRARIES=/ucrt64/lib/libmsmpi.a \
  -DMPI_GUESS_LIBRARY_NAME=MSMPI .

mingw32-make