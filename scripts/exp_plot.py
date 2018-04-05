import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib2tikz import save as tikz_save
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import re
from scipy.stats import norm, stats
import uncertainties as unc
import os
import glob
import re
import linecache
mpl.style.use('seaborn')
sns.set_style("white")
mpl.rcParams.update({'font.size': 32})

regU = r"U([0-9]+\.[0-9]+)"
regb = r"b\_([0-9]+\.[0-9]+)"
regSeMe = r"iExpansionOrder(.+)\.out$"
path = "./ExpansionOrder*"
data = []

for filename in glob.glob(path):
    with open(filename, 'r') as infile:
        U = float(re.search(regU, filename).group(0)[1:])
        b = float(re.search(regb, filename).group(0)[2:])
        _ = infile.readline()
        _ = infile.readline()
        tmp = np.fromstring(infile.readline(), sep=' ')
        data.append([b, U, 0.25 - tmp[1]/(b*U) + (0.502*0.502), 1. - 2* (0.5 - tmp[1]/(b*U)), np.sqrt(tmp[2]/tmp[0])])
data = np.array(data)
d10 = data[data[:,0] == 10]
d10s = d10[d10[:,1].argsort()]
d15 = data[data[:,0] == 15]
d15s = d15[d15[:,1].argsort()]
d20 = data[data[:,0] == 20]
d20s = d20[d20[:,1].argsort()]
d25 = data[data[:,0] == 25]
d25s = d25[d25[:,1].argsort()]
d30 = data[data[:,0] == 30]
d30s = d30[d30[:,1].argsort()]
plt.errorbar(d10s[:,1], d10s[:,2], yerr=d10s[:,4],  label=r"$\beta = 10$")
plt.errorbar(d15s[:,1], d15s[:,2], yerr=d15s[:,4], label=r"$\beta = 15$")
plt.errorbar(d20s[:,1], d20s[:,2], yerr=d20s[:,4], label=r"$\beta = 20$")
plt.errorbar(d25s[:,1], d25s[:,2], yerr=d25s[:,4],  label=r"$\beta = 25$")
plt.errorbar(d30s[:,1], d30s[:,2], yerr=d30s[:,4],  label=r"$\beta = 30$")
plt.axvline(x=2.3,ymax = (d10s[5,2]/plt.ylim()[1]), ls = '-.', alpha = 0.7)
plt.ylabel(r"$\langle n_\uparrow n_\downarrow \rangle$")
plt.xticks(list(plt.xticks()[0]) + [2.3])
plt.xlabel("U")
plt.legend()
plt.show()
plt.plot(d10s[:,1], d10s[:,3], '-o', label=r"$\beta = 10$")
plt.plot(d15s[:,1], d15s[:,3], '-o', label=r"$\beta = 15$")
plt.plot(d20s[:,1], d20s[:,3], '-o', label=r"$\beta = 20$")
plt.plot(d25s[:,1], d25s[:,3], '-o', label=r"$\beta = 25$")
plt.plot(d30s[:,1], d30s[:,3], '-o', label=r"$\beta = 30$")
plt.ylabel(r"$\Delta n$")
plt.xlabel("U")
plt.legend()
plt.show()
plt.semilogy(d10s[:,1], d10s[:,3], '-o', label=r"$\beta = 10$")
plt.semilogy(d15s[:,1], d15s[:,3], '-o', label=r"$\beta = 15$")
plt.semilogy(d20s[:,1], d20s[:,3], '-o', label=r"$\beta = 20$")
plt.semilogy(d25s[:,1], d25s[:,3], '-o', label=r"$\beta = 25$")
plt.semilogy(d30s[:,1], d30s[:,3], '-o', label=r"$\beta = 30$")
plt.xlabel("U")
plt.legend()
plt.show()
