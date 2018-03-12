#!/usr/bin/env python
import os.path
import glob
import re
import numpy as np
from sys import argv
from time import clock
from subprocess import call

regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
regSeMe = r"f_SelfE(.+)[0-9]+\.out$"
regGImpMe = r"f_GImp(.+)[0-9]+\.out$"

#from: https://stackoverflow.com/questions/1557571/how-do-i-get-time-of-a-python-programs-execution
def timeFormat(t):
	return "%d:%02d:%02d.%03d" % reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t*1000,),1000,60,60])

def getOpts(argv):
	opts = {}
	while argv:
		if argv[0][0] == '-':
			if len(argv) > 1:
				opts[argv[0]] = argv[1]
			else:
				opts[argv[0]] = True
		argv = argv[1:]
	return opts

start_time = clock()
i = 0
tot =0
silent = (True if "-silent" in getOpts(argv) else False)
with open(os.devnull, 'w') as devnull:
	for filename in glob.glob("me_GImp*.in"):
	    tot += 1
	    if not os.path.isfile(filename[:-2] + "out.avspec.dat"):
		i += 1

	ii = 0
	for filename in glob.glob("me_GImp*.in"):
		if os.path.isfile(filename[:-2] + "out.avspec.dat"):
			pass		
			#print("Skipping " + filename + " because " + filename[:-2] + "out.avspec.dat" + " exists")
		else:
			ii += 1
			print("[" + timeFormat(clock()-start_time) + "] running file " + str(ii) + "/" +str(i) + " from a total of " +str(tot))
			if silent:
				call(["maxent", filename], stdout=devnull, stderr=devnull)
			else:
				print("Filename: " + filename)
				call(["maxent", filename])

