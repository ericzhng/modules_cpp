
import os
dllLoc = os.path.abspath("./build/Release")

import sys
sys.path.append(dllLoc)
print(sys.path)

import myboostmodule

g = myboostmodule.Greeter()
g.say_hello("Bob")
