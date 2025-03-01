
import ctypes
from cffi import FFI

import struct
print(struct.calcsize("P") * 8)  # Outputs 32 or 64

import os
# dllLoc = os.path.abspath("./build/Debug/mycpp.dll")
dllLoc = os.path.abspath("./build/Release/mycpp.dll")



# Create a FFI instance
ffi = FFI()
# Load the DLL (replace 'mydll.dll' with the actual path to your DLL)
try:
    mycpp = ffi.dlopen(dllLoc)
except OSError as e:
    print(f"Error loading DLL: {e}")
    exit()

ffi.cdef("""
    void hello(int x);
    """)

# Call the C function
mycpp.hello(42)




# Load the shared library
mycpp1 = ctypes.CDLL(dllLoc)

print(dir(mycpp1))

# Declare function argument and return types
mycpp1.hello.argtypes = [ctypes.c_int]
mycpp1.hello.restype = None  # Function returns void

# Call the function
mycpp1.hello(ctypes.c_int(42))




mycpp2 = ctypes.windll.LoadLibrary(dllLoc)

# Access a function from the DLL
mycpp2.hello(ctypes.c_int(42))
