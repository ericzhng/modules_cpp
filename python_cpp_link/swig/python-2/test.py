
import os
dllLoc = os.path.abspath("./build/Release")
print(dllLoc)

import sys
sys.path.append(dllLoc)

import mycpp2

greeter = mycpp2.Greeter()
greeter.say_hello("Alice")
