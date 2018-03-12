#!/usr/bin/env python
import os.path
import glob
import re
import numpy as np
from subprocess import call

regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
regSeMe = r"f_SelfE(.+)[0-9]+\.out$"
regGImpMe = r"f_GImp(.+)[0-9]+\.out$"
i = 0
tot =0
for filename in glob.glob("me_GImp*.in"):
    tot += 1
    if os.path.isfile(filename[:-2] + "out.avspec.dat"):
        pass
        #print("Skipping " + filename + " because " + filename[:-2] + "out.avspec.dat" + " exists")
    else:
        i += 1
        #print("running: " + filename)
        #call(["maxent", filename])
print(str(i) + " of " + str(tot) +" files not continued")
