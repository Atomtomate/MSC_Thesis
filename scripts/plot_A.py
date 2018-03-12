import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
from scipy.stats import norm
import os
import glob
import re
import linecache

mpl.style.use('seaborn')
sns.set_style("whitegrid")
regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
data = []

file_list = ["me_GImp_b100_000000_U2_650000conf.out.avspec.dat", "me_GImp_b100_000000_U2_900000conf.out.avspec.dat","me_GImp_b100_000000_U3_400000conf.out.avspec.dat"]
for filename in file_list:
    with open(filename, 'r') as infile:
	d = np.genfromtxt(infile)
	U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
	ind = np.where(np.abs(d[:,0] < 1.5*U))
	#dd = d[ind]
	data.append(d)
	plt.plot(d[:,0], d[:,1], label="U/D = "+str(U))
	plt.fill_between(d[:,0], 0, d[:,1], alpha=0.3)
#plt.ylim((0.,))
plt.xlim([-3.5,3.5])
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")
plt.legend()
plt.show()
