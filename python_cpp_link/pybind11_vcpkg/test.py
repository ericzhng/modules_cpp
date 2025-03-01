
import os
dllLoc = os.path.abspath("./build/Release")
print(dllLoc)

import sys
sys.path.append(dllLoc)


import myadd
print(myadd.add(21, 3))  # Output: 5


import myadd2
print(myadd2.add(21, 3))  # Output: 5
