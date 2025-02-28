
import os
import ctypes

dllLoc = os.path.abspath("./build/Debug/ctypes.dll")

# Load the shared library
mycpp = ctypes.CDLL(dllLoc)

print(dir(mycpp))

# Declare function argument and return types
mycpp.hello.argtypes = [ctypes.c_int]
mycpp.hello.restype = None  # Function returns void

# Call the function
mycpp.hello(ctypes.c_int(42))


mycpp2 = ctypes.windll.LoadLibrary(dllLoc)
# Access a function from the DLL
mycpp2.hello(ctypes.c_int(42))
