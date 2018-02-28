import numpy as np
import os
import glob

for filename in glob.glob('*out.avspec.dat'):
	with open(filename, 'r') as infile:
		arr = np.genfromtxt(infile)
		print(filename)
		print(arr)

