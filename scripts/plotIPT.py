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
regU = r"U([0-9]+\.[0-9]+)"
data = []
for filename in glob.glob(r'*pade.cont'):
    with open(filename, 'r') as infile:
	d = np.genfromtxt(infile,skip_footer=1)
	U = float(re.search(regU, filename).group(0)[1:])
	ind = np.where(np.abs(d[:,0] < 1.5*U))
	#dd = d[ind]
	data.append(d)
	if(U < 6):
		plt.plot(d[:,0], -d[:,2]/np.pi, label="U = "+str(U))
		plt.fill_between(d[:,0], 0, -d[:,2]/np.pi, alpha=0.3)
plt.ylim((0.,0.7))
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")
plt.legend()
plt.show()
