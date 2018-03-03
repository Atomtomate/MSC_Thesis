#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.stats import norm
import os
import glob
import re
import linecache

mpl.style.use('seaborn')
sns.set_style("whitegrid")

def plot_exp(data, title):
    xticks = [0., 3., 22.]
    for i in range(len(data)):
        h = data_hist[i]
        d = data[i]
        if(data[i][3] > 0.8 or data[i][4] > 0.8):
            print("Not a normal distribution!")
            print(data[i])
        x = np.linspace(0, data[i][1]+1.2*data[i][2], 3000)
        if(d[0]%1 == 0):
            p = plt.plot(x, norm.pdf(x, data[i][1], np.sqrt(data[i][2])) , label="U=" + str(d[0]))
            ax = plt.gca()
            ax.set_ylim(bottom=0.,top=0.25)
            ax.axvline(data[i][1], ymax=(norm.pdf(data[i][1], data[i][1], np.sqrt(data[i][2]))/(ax.get_ylim()[1])), alpha=0.4, c=p[-1].get_color())
            xticks.append(data[i][1])
    plt.ylabel(r"$p(\langle n \rangle )$")
    plt.xlabel(r"$\langle n \rangle$")
    plt.title(title)
    xticks_l = [ '%.1f' % el for el in xticks ]
    plt.xticks(xticks, xticks_l)
    plt.grid(False)
    plt.legend()
    plt.show()


dataH = []
dataH_hist = []
for filename in glob.glob(r'./expOrder_CTHYB/*ExpansionOrder*'):
    with open(filename, 'r') as infile:
        print(filename)
        t = re.search("U([0-9]+\.[0-9]+)\.out", filename)
        tmp = [float(t.group(0)[1:-4])]
        tmp.extend(map(float, linecache.getline(filename, 2).split()))
        tmp_hist = np.genfromtxt(infile, skip_header = 3)
        tmp_hist = tmp_hist[tmp_hist[:,1] != 0]
        dataH.append(tmp)
        dataH_hist.append(tmp_hist)
dataH.sort(key=lambda x: x[0])
plt.figure(0)
#plot_exp(dataH, "Expansion Order for CT-HYB")

dataI = []
dataI_hist = []
for filename in glob.glob(r'./expOrder_CTINT/*ExpansionOrder*'):
    with open(filename, 'r') as infile:
        print(filename)
        t = re.search("U([0-9]+\.[0-9]+)\.out", filename)
        tmp = [float(t.group(0)[1:-4])]
        tmp.extend(map(float, linecache.getline(filename, 2).split()))
        tmp_hist = np.genfromtxt(infile, skip_header = 3)
        tmp_hist = tmp_hist[tmp_hist[:,1] != 0]
        dataI.append(tmp)
        dataI_hist.append(tmp_hist)
dataI.sort(key=lambda x: x[0])
plt.figure(1)
#plot_exp(dataI, "Expansion Order for CT-INT")

plt.figure(2)
yH = []
yeH = []
xH = []
yI = []
yeI = []
xI = []
for el in dataH:
    xH.append(el[0])
    yH.append(el[1])
    yeH.append(np.sqrt(el[2]))
for el in dataI:
    xI.append(el[0])
    yI.append(el[1])
    yeI.append(np.sqrt(el[2]))
print(xH)
print(yH)
plt.title(r"Expansion order at $\beta = 32$")
plt.errorbar(xH, yH, yerr=yeH, linestyle='-', marker='o', label="CT-HYB")
plt.errorbar(xI, yI, yerr=yeI, linestyle='-', marker='o', label="CT-INT")
plt.xlabel("U")
plt.ylabel(r"$\langle n \rangle$")
plt.legend()
plt.show()
