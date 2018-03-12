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

for filename in glob.glob("me_GImp*.in"):
	if os.path.isfile(filename[:-2] + "out.avspec.dat"):
		print("Skipping " + filename + "\nbecause " + filename[:-2] + "out.avspec.dat" + "exists")
	else:
		print("running: " + filename)
		call(["maxent", filename])

