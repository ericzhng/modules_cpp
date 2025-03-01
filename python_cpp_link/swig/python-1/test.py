
import os
dllLoc = os.path.abspath("./build/Release")
print(dllLoc)

import sys
sys.path.append(dllLoc)

import mycpp

greeter = mycpp.Greeter()
greeter.say_hello("Alice")
