#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
from scipy.stats import norm, stats
import os
import glob
import re
import linecache

mpl.style.use('seaborn')
sns.set_style("whitegrid")

use_pade = False
use_me = True
regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
regSeMe = r"f_SelfE(.+)[0-9]+\.out$"
data = []
se_data = []
pade_str = r'*pade.cont'
me_str = r'*avspec.dat'
se_str = r'f_SelfE*0.out'

""" returns the width of suppress spectral weight around the fermi level
"""
def comp_d(data, epsilon=0.01):
    mi = np.argmax(dat[:,0]>0)
    i = 1
    d = 0.
    sweight = 0.
    while mi - i >= 0:
        dat_n = dat[mi-i:mi+i,:]
        d_n = dat_n[-1,0] - dat_n[0,0]
        sweight_n = np.sum(dat_n[:,1])
        if sweight_n > epsilon:
            break
        d = d_n
        sweight = sweight_n
        i += 1
    return d

def comp_Z(se_data):
    ulist = np.unique(se_data[:,1])
    max_points = 3
    Z = []
    for u in ulist:
        ui =  np.where(se_data[:,1] == u)
        #find lowest available temperatures
	ii = np.argsort(se_data[ui][:,0])
        d = se_data[ui][ii][-max_points:]	# list of lowest temperatur for given U
	w0l = np.pi/d[:,0]			# zero Matsubara Frequency
	dRSigma = d[:,2]/w0l			# approximation for the derivative of SE at w=0
	dRSigma = stats.linregress(1./np.array(d[:,0]), dRSigma).intercept
	Z.append([u,1./(1.-dRSigma)])
    return np.array(Z)

""" compute critical value of U for all temparatures.

    data is a array containing arrays with 3 entries: [beta, U, d].
    d is obtained by calling comp_d
"""
def comp_pt(data):
    beta_list = np.unique(data[:,0])
    for beta in beta_list:
        d20 = np.array(data[np.where(data == beta)[0]])
        d20s = np.array(sorted(d20, key=lambda el : el[1]))


if use_pade:
    for filename in glob.glob(pade_str):
        with open(filename, 'r') as infile:
            d = np.genfromtxt(infile)
            dd = d[:-1,:]
            info = d[-1,2]
            insulating = 1
            if(info < -0.1):
                    insulating = -1
            U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
            b = float(re.search(regb, filename).group(0)[1:].replace("_","."))
            data.append([b, U, insulating])
if use_me:
    for filename in glob.glob(me_str):
        with open(filename, 'r') as infile:
            dat = np.genfromtxt(infile)
            U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
            b = float(re.search(regb, filename).group(0)[1:].replace("_","."))
            d = comp_d(dat)
            data.append([b,U,d])
print("done reading MaxEnt data\n")
i = 0
for filename in glob.glob(se_str):
    with open(filename, 'r') as infile:
	infile.readline()
        dat = np.fromstring(infile.readline(), sep=' ')
	U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
        b = float(re.search(regb, filename).group(0)[1:].replace("_","."))
        se_data.append([b,U,dat[3]])
print("done reading Self Energy data\n")
data = np.array(data)
se_data = np.array(se_data)
Zlist = comp_Z(se_data)

nb = np.unique(data[:,0]).shape[0]
bSorted = np.unique(data[:,0])
bSorted.sort()
uSorted = np.unique(data[:,1])
uSorted.sort()
nU = np.unique(data[:,1]).shape[0]
image = 255*np.ones((nU,nb))
for el in data:
	yi = np.where(bSorted == el[0])[0][0]
	xi = np.where(uSorted == el[1])[0][0]
	image[xi, -yi] = el[2]

plt.imshow(image.T, interpolation='bilinear', cmap='viridis')
ax = plt.gca()
ax.set_xticklabels(uSorted)
ax.set_yticklabels(1./bSorted)
plt.show()

