
import os
dllLoc = os.path.abspath("../build/Release")
print(dllLoc)

import sys
sys.path.append(dllLoc)

import swig_ex
swig_ex.say_hello()
