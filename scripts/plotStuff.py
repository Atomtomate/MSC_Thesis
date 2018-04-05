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
mpl.rcParams.update({'font.size': 32})
mpl.style.use('seaborn')
sns.set_style("whitegrid")
regU = r"U([0-9]+\_[0-9]+)"
data = []
dataSE = []
files = ["me_GImp_b30_000000_U2_200000conf.out.avspec.dat", "me_GImp_b30_000000_U2_300000conf.out.avspec.dat", "me_GImp_b30_000000_U2_400000conf.out.avspec.dat", "me_GImp_b30_000000_U2_500000conf.out.avspec.dat"]
filesSE = ["f_SelfE_b30_000000_U2_200000.out", "f_SelfE_b30_000000_U2_300000.out", "f_SelfE_b30_000000_U2_400000.out", "f_SelfE_b30_000000_U2_500000.out"]
filesGI = ["f_GImp_b30_000000_U2_200000.out", "f_GImp_b30_000000_U2_300000.out", "f_GImp_b30_000000_U2_400000.out", "f_GImp_b30_000000_U2_500000.out"]
for filename in files:#glob.glob(r'*pade.cont'):
    with open(filename, 'r') as infile:
	d = np.genfromtxt(infile,skip_footer=1)
	U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
	ind = np.where(np.abs(d[:,0] < 1.5*U))
	#dd = d[ind]
	data.append(d)
	plt.plot(d[:,0], d[:,1], label="U = "+str(U))
	plt.fill_between(d[:,0], 0, d[:,1], alpha=0.3)
	#plt.plot(d[:,0], -d[:,2]/np.pi, label="U = "+str(U))
	#plt.fill_between(d[:,0], 0, -d[:,2]/np.pi, alpha=0.3)
#plt.ylim((0.,0.7))
plt.xlim((-3,3))
plt.xlabel(r"$\omega$")
plt.ylabel(r"$A(\omega)$")
plt.legend()
plt.show()
for filename in filesSE:
    with open(filename, 'r') as infile:
	d = np.genfromtxt(infile,skip_header=1)
	U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
	ind = np.where(np.abs(d[:,0] < 1.5*U))
	#dd = d[ind]
	data.append(d)
        plt.errorbar(d[:20,0], d[:20,3], fmt='.' , yerr=d[:20,4] , label="U = "+str(U))
#plt.ylim((0.,0.7))
plt.xlabel(r"$i \omega_n$")
plt.ylabel(r"$\Sigma (i \omega_n)$")
#plt.ylabel(r"$G_{imp}(i \omega_n)$")
plt.legend()
plt.show()
