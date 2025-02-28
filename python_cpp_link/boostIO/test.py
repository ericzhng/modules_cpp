
import os
dllLoc = os.path.abspath("./build/Release")

import sys
sys.path.append(dllLoc)
import sys
print(sys.path)

import boost_python

g = boost_python.Greeter()
g.say_hello("Bob")
