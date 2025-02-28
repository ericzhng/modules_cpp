
import struct
print(struct.calcsize("P") * 8)  # Outputs 32 or 64

import os
dllLoc = os.path.abspath("./build/Debug/ffi.dll")

from cffi import FFI
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
