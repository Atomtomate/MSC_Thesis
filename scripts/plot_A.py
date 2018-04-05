import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import seaborn as sns
import re
from scipy.stats import norm
import os
import glob
import re
import linecache

mpl.style.use('seaborn')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
sns.set_style("whitegrid")
regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
data = []

file_list =["me_GImp_b50_000000_U2_390000conf.out.avspec.dat", "me_GImp_b50_000000_U2_400000conf.out.avspec.dat", "me_GImp_b50_000000_U2_410000conf.out.avspec.dat", "me_GImp_b50_000000_U2_420000conf.out.avspec.dat", "me_GImp_b50_000000_U2_500000conf.out.avspec.dat", "me_GImp_b50_000000_U3_200000conf.out.avspec.dat"]
plot_x = []
plot_y = []
for filename in file_list:
    with open(filename, 'r') as infile:
	d = np.genfromtxt(infile)
	U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
	ind = np.where(np.abs(d[:,0] < 1.5*U))
	#dd = d[ind]
	data.append(d)
        plot_x.append(d[:,0])
        plot_y.append(d[:,1])
        plt.plot(d[:,0], d[:,1], label="U/D = "+str(U))
        plt.fill_between(d[:,0], 0., d[:,1] , alpha=0.3)
plot_x = np.array(plot_x).T
plot_y = np.array(plot_y).T
z_x = np.zeros(plot_x.shape[1])

#plt.ylim((0.,))
plt.xlim([-3.0,3.0])
plt.xlabel(r"$\omega$", fontsize=18)
plt.ylabel(r"$A(\omega)$", fontsize=18)
plt.legend(fontsize=18)
ax = plt.gca()
axins = inset_axes(ax,1.8,1.8, loc=1, bbox_to_anchor=(0.28, 0.85), bbox_transform=ax.figure.transFigure)
axins.set_xlim(-0.23,0.23)
axins.set_ylim(0.0,0.0035)
plt.plot(plot_x, plot_y, label="U/D = "+str(U))
box, c1, c2 = mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
plt.setp([c1, c2], color="#4a4a4a", linewidth=1.2, linestyle=":")
plt.setp(box, linewidth=1, color="#4a4a4a", linestyle=":")
plt.show()
