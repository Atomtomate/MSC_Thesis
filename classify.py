#!/usr/bin/env python

import matplotlib as mpl
#mpl.use('pgf')
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
cmap = mpl.cm.GnBu

use_pade = False
use_me = True
regU = r"U([0-9]+\_[0-9]+)"
regb = r"b([0-9]+\_[0-9]+)"
regSeMe = r"f_SelfE(.+)[0-9]+\.out$"
data = []
pade_str = r'*pade.cont'
me_FROM0_str = r'IPT_Bethe_PD_from_0/*avspec.dat'
me_TO0_str = r'IPT_Bethe_PD_to_0/*avspec.dat'
se04_str = r'IPT_Bethe_PD_from_0/f_SelfE*0.out'
se40_str = r'IPT_Bethe_PD_to_0/f_SelfE*0.out'
se23_str = r'IPT_PD_Z_03/f_SelfE*0.out'
se32_str = r'IPT_PD_Z_30/f_SelfE*0.out'

""" returns the width of suppress spectral weight around the fermi level
"""
def comp_d(data, epsilon=1e-15, epsilon2=0.1):
    mi = np.argmax(data[:,0]>0)
    i = 1
    d = 0.
    sweight = 0.
    while mi - i >= 0:
        dat_n = data[mi-i:mi+i,:]
        d_n = dat_n[-1,0] - dat_n[0,0]
        sweight_n = np.sum(dat_n[:,1])
        if sweight_n > epsilon:
            break
        d = d_n
        sweight = sweight_n
        i += 1
    if d < epsilon2:
	return 0.
    return d

def comp_Z0(se_data):
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
	res = stats.linregress(1./np.array(d[:,0]), dRSigma)
	rr = unc.ufloat(res.intercept,res.stderr)
	rr = (1./(1.-rr))
	Z.append([u,rr.n,rr.std_dev])
    return np.array(Z)

def comp_Z(se_data):
    dRSigma = se_data[:,2]
    w0l = np.pi/se_data[:,0]
    rr = (1./(1.-dRSigma/w0l))
    return np.hstack((se_data, rr[:,np.newaxis]))

""" compute critical value of U for all temparatures.

    data is a array containing arrays with 3 entries: [beta, U, d].
    d is obtained by calling comp_d
"""
def comp_pt(data):
    beta_list = np.unique(data[:,0])
    for beta in beta_list:
        d20 = np.array(data[np.where(data == beta)[0]])
        d20s = np.array(sorted(d20, key=lambda el : el[1]))

def read_se(path, nPoints=1):
    se_data = []
    for filename in glob.glob(path):
        with open(filename, 'r') as infile:
	    U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
            b = float(re.search(regb, filename).group(0)[1:].replace("_","."))
	    last_pos = infile.tell()
	    line = infile.readline()
	    if not (line.split()[0] == "mFreq"):
		infile.seek(last_pos)
	    dat = []
	    for i in range(nPoints):
                dat.append(np.fromstring(infile.readline(), sep=' ')[3])
            se_data.append([b,U]+dat)
    return np.array(se_data)

def read_A(path):
    data = []
    if use_pade:
        for filename in glob.glob(path):
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
        for filename in glob.glob(path):
            with open(filename, 'r') as infile:
                dat = np.genfromtxt(infile)
                U = float(re.search(regU, filename).group(0)[1:].replace("_","."))
                b = float(re.search(regb, filename).group(0)[1:].replace("_","."))
                d = comp_d(dat)
                data.append([b,U,d])
    return np.array(data)


se_data_to = read_se(se32_str,1)#se40_str, 1)

#data_from = read_A(me_FROM0_str)
print("read from MI to metal")
#data_to = read_A(me_TO0_str)
se_data_from = read_se(se23_str,)#se04_str, 1)
print("read from metal to MI")



def plotZ():
    se04_data = read_se(se23_str)#se04_str)
    se40_data = read_se(se32_str)#se40_str)
    Z04list = comp_Z0(se04_data)
    Z40list = comp_Z0(se40_data)
    y04 = np.clip(Z04list[:,1],a_min=0.00001, a_max=None)
    y40 = np.clip(Z40list[:,1],a_min=0.00001, a_max=None)
    #pti = np.argmax(y < 0.0001)
    #yerr = np.clip(Zlist[:,2], a_min=0.,a_max=1.)
    #yerr[pti:] = 0.
    #y[pti:] = 0.
    plt.semilogy(Z40list[:,0], y40, "o-", ms=3.0, markevery=4, label="init U = 4")
    plt.semilogy(Z04list[:,0], y04, 'o-', alpha=0.4, ms=4.0, markevery=4, label="init U = 0")
    plt.xlabel(r"U")
    plt.ylabel(r"Z")
    plt.legend()
    plt.show() 
    plt.plot(Z40list[:,0], y40, "o", ms=3.0, label="init U = 4")
    plt.plot(Z04list[:,0], y04, 'o', alpha=0.4, ms=4.0, label="init U = 0")
    plt.xlabel(r"U")
    plt.ylabel(r"Z")
    plt.legend()
    plt.show()


def tmp():
	data = data_to
	blist = np.unique(data[:,0])
	ptl = []
	for beta in blist:
	    bi = np.where(data[:, 0] == beta)
	    dat = np.squeeze(data[bi,:])
	    dat = dat[dat[:,1].argsort()]
	    pti = (dat[:,2] != 0.).argmax()
	    ptl.append([beta, dat[pti,1]])
	ptl = np.array(ptl)
	invalid = 4#(ptl[:, 1] == 2.).argmin()
	fig, ax = plt.subplots()
	ax.set_xlim([2.2, 4.0])
	ax.xaxis.set_major_locator(MaxNLocator(8))
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
	ax.plot(ptl[invalid:,1], 1./ptl[invalid:,0], label="to")

	data = data_from

	blist = np.unique(data[:,0])
	ptl = []
	for beta in blist:
	    bi = np.where(data[:, 0] == beta)
	    dat = np.squeeze(data[bi,:])
	    dat = dat[dat[:,1].argsort()]
	    pti = (dat[:,2] != 0.).argmax()
	    ptl.append([beta, dat[pti,1]])
	ptl = np.array(ptl)
	invalid = 4#(ptl[:, 1] == 2.).argmin()


	ax.plot(ptl[invalid:,1], 1./ptl[invalid:,0], label="from")
	plt.show()

""" returns [beta, U, Im(0), d Im(eps)/ d w, fermi_liquid_check]
    fermi_liquid_check == 1 if Im(0) < eps2 && d Im(eps) / d w < eps2 else 0.
    The values of the imeginary part for the self energy are obtained by a
    polynomial fit of order 4 to the lowest nPoints Matsubara frequencies
"""
def checkLF(nPoints= 6, eps = 0.001, eps2 = 0.05):
    seT = read_se(se40_str, 6)
    seF = read_se(se04_str, 6)
    res_BU = []
    for data in [seT, seF]:
	nb = np.unique(data[:,0]).shape[0]
	nU = np.unique(data[:,1]).shape[0]
	blist = np.unique(data[:,0])
	ulist = np.unique(data[:,1])
	blist.sort()
	ulist.sort()
	res_line = np.zeros((nb, nU, 3))
	for el in data:
		wnlist = np.array([1j*(2.*n+1.)*np.pi/el[0] for n in range(len(el)-2)])
		res = np.polyfit(wnlist, 1j*el[2:], 4).imag
		res = np.poly1d(res)
		resd = np.polyder(res)
		fl_check = (1. if np.polyval(res,0.) < eps2 and  np.polyval(resd,eps) < eps2 else 0.)
		bi = np.squeeze(np.where(blist == el[0]))
		ui = np.squeeze(np.where(ulist == el[1]))
		res_line[bi,ui,:] = np.array([np.polyval(res,0.), np.polyval(resd,eps), fl_check], dtype=np.float64)
		#res_line.append([el[0], el[1], res(0.), res(eps), fl_check])
        res_BU.append(res_line)
    res_BU = np.array(res_BU)
    fig, axarr = plt.subplots(2,2)
    fig.subplots_adjust(wspace=.4, hspace=0.3)
    d = np.clip(res_BU[0,:,:,0],None, 0.) #clip unphysical negative values, they are noise from simulation
    vmax = np.max(d) + 0.01
    cax1 = axarr[0,0].imshow(-d+vmax, norm=colors.LogNorm(vmin=-np.min(d)+vmax, vmax=-np.max(d)+vmax), interpolation='None', cmap=cmap, aspect='auto')
    axarr[0,0].set_title(r"$\log(-Im(\Sigma(0)))$, IG $U=4$")

    vmin = np.min(res_BU[0,:,:,0]) + 0.01
    cax2 = axarr[0,1].imshow(res_BU[0,:,:,1]-vmin, norm=colors.LogNorm(vmin=np.min(res_BU[0,:,:,1])-vmin, vmax=np.max(res_BU[0,:,:,1])-vmin), interpolation='None', cmap=cmap, aspect='auto')
    axarr[0,1].set_title(r"$\frac{\partial\, Im[\Sigma(\omega)]}{\partial \omega}$, IG $U=4$")

    d = np.clip(res_BU[1,:,:,0],None, 0.) #clip unphysical negative values, they are noise from simulation
    vmax = np.max(d) + 0.01
    cax3 = axarr[1,0].imshow(-d+vmax, norm=colors.LogNorm(vmin=-np.min(d)+vmax, vmax=-np.max(d)+vmax), interpolation='None', cmap=cmap, aspect='auto')
    axarr[1,0].set_title(r"$\log(-Im(\Sigma(0)))$, IG $U=0$")

    vmin = np.min(res_BU[1,:,:,0]) + 0.01
    cax4 = axarr[1,1].imshow(res_BU[1,:,:,1]-vmin, norm=colors.LogNorm(vmin=np.min(res_BU[1,:,:,1])-vmin, vmax=np.max(res_BU[1,:,:,1])-vmin), interpolation='None', cmap=cmap, aspect='auto')
    axarr[1,1].set_title(r"$\frac{\partial\, Im[\Sigma(\omega)]}{\partial \omega}$, IG $U=0$")
    for i in range(axarr.shape[0]):
	    for ax in axarr[i,:]:
		ax.set_yticks(np.arange(0,len(blist),5))
		ax.set_yticklabels(blist[::5].astype(np.int))
		ax.set_xticks(np.arange(0,len(ulist),6))
		ax.set_xticklabels(ulist[::6])
		ax.set_xlabel("U")
		ax.set_ylabel(r"$\beta$")
    divider = make_axes_locatable(axarr[0,0])
    cax = divider.append_axes('right', size=0.2, pad=0.05)
    fig.colorbar(cax1, cax=cax, orientation='vertical')
    divider = make_axes_locatable(axarr[0,1])
    cax = divider.append_axes('right', size=0.2, pad=0.05)
    fig.colorbar(cax2, cax=cax, orientation='vertical')
    divider = make_axes_locatable(axarr[1,0])
    cax = divider.append_axes('right', size=0.2, pad=0.05)
    fig.colorbar(cax3, cax=cax, orientation='vertical')
    divider = make_axes_locatable(axarr[1,1])
    cax = divider.append_axes('right', size=0.2, pad=0.05)
    fig.colorbar(cax4, cax=cax, orientation='vertical')
    tikz_save("FL_check1.pgf")
    fig,axarr = plt.subplots(1,2)
    axarr[0].imshow(-res_BU[0,:,:,2], cmap=mpl.cm.binary, interpolation='spline36', aspect='auto')
    axarr[0].set_title(r"Valid LF, initial guess $U=4$")
    axarr[0].set_yticks(np.arange(0,len(blist),5))
    axarr[0].set_yticklabels(blist[::5].astype(np.int))
    axarr[0].set_xticks(np.arange(0,len(ulist),6))
    axarr[0].set_xticklabels(ulist[::6])
    axarr[0].set_xlabel("U")
    axarr[0].set_ylabel(r"$\beta$")
    axarr[1].imshow(-res_BU[1,:,:,2], cmap=mpl.cm.binary, interpolation='spline36', aspect='auto')
    axarr[1].set_title(r"Valid LF, initial guess $U=0$")
    axarr[1].set_yticks(np.arange(0,len(blist),5))
    axarr[1].set_yticklabels(blist[::5].astype(np.int))
    axarr[1].set_xticks(np.arange(0,len(ulist),6))
    axarr[1].set_xticklabels(ulist[::6])
    axarr[1].set_xlabel("U")
    axarr[1].set_ylabel(r"$\beta$")
    plt.show()
    tikz_save("FL_check2.pgf")
    #plt.savefig("FL_check2.pgf")
    return res_BU
		

def plotPD():
	data_to0=se_data_to
	data_from0=se_data_from
	val_axis=3
	#def plotA(data_to0, data_from0, val_axis=2):
	dt0 = comp_Z(data_to0)
	df0 = comp_Z(data_from0)
	#cmap.set_bad(color='black')
	nb = np.unique(dt0[:,0]).shape[0]
	nU = np.unique(dt0[:,1]).shape[0]
	blist = np.unique(dt0[:,0])
	ulist = np.unique(dt0[:,1])
	blist.sort()
	ulist.sort()
	imt0 = np.ones((nU,nb))
	for bi,beta in enumerate(blist):
		bii = np.where(dt0[:,0] == beta)
		for ui,U in enumerate(ulist):
                    uii = np.where(dt0[:,1] == U)
                    di = np.intersect1d(bii,uii)
                    imt0[ui,bi] = dt0[di, val_axis]
	# maybe don't copy paste so much?...
	nb = np.unique(df0[:,0]).shape[0]
	nU = np.unique(df0[:,1]).shape[0]
	blist = np.unique(df0[:,0])
	ulist = np.unique(df0[:,1])
	blist.sort()
	ulist.sort()
	imf0 = np.ones((nU,nb))
	for bi,beta in enumerate(blist):
		bii = np.where(df0[:,0] == beta)
		for ui,U in enumerate(ulist):
                    uii = np.where(df0[:,1] == U)
                    di = np.intersect1d(bii,uii)
                    imf0[ui,bi] = df0[di, val_axis]
	imdiff = imf0.T-imt0.T
	imt0 = imt0.T
	imf0 = imf0.T
	#imdiff = np.ma.masked_where(imdiff < 0.01, imdiff)
	#imt0 = np.ma.masked_where(imt0.T < 0.01, imt0.T)
	#imf0 = np.ma.masked_where(imf0.T < 0.01, imf0.T)
	 

	fig,ax = plt.subplots(1,2)
	iax = ax[0].imshow(imt0, interpolation='nearest', cmap=cmap, aspect='auto')
	divider = make_axes_locatable(ax[1])
 	cax = divider.append_axes('right', size=0.2, pad=0.05)
	fig.colorbar(iax, cax=cax, orientation='vertical')
	ax[0].set_yticks(np.arange(0,len(blist),3))
	ax[0].set_yticklabels(blist[::3].astype(np.int))
	ax[0].set_xticks(np.arange(0,len(ulist),3))
	ax[0].set_xticklabels(ulist[::3])
	ax[0].set_xlabel("U")
	ax[0].set_ylabel(r"$\beta$")
	ax[0].set_title(r"$Z(\beta, U)$, initial guess $U=0$")

	#fig,ax = plt.subplots()
	iax = ax[1].imshow(imf0, interpolation='nearest', cmap=cmap, aspect='auto')
	#, extent=[ulist[0], ulist[-1], blist[-1], blist[0]]
	ax[1].set_yticks(np.arange(0,len(blist),3))
	ax[1].set_yticklabels(blist[::3].astype(np.int))
	ax[1].set_xticks(np.arange(0,len(ulist),3))
	ax[1].set_xticklabels(ulist[::3])
	divider = make_axes_locatable(ax[1])
 	cax = divider.append_axes('right', size=0.2, pad=0.05)
	fig.colorbar(iax, cax=cax, orientation='vertical')
	ax[1].set_xlabel("U")
	ax[1].set_ylabel(r"$\beta$")
	ax[1].set_title(r"$Z(\beta, U)$, initial guess $U=4$")

	fig,ax = plt.subplots()
	cax = ax.imshow(imdiff, interpolation='nearest', cmap=cmap)
	#, extent=[ulist[0], ulist[-1], blist[0], blist[-1]]
	ax.set_yticks(np.arange(0,len(blist),3))
	ax.set_yticklabels(blist[::3].astype(np.int))
	ax.set_xticks(np.arange(0,len(ulist),3))
	ax.set_xticklabels(ulist[::3])
	#divider = make_axes_locatable(ax)
 	#iax = divider.append_axes('right', size=0.2, pad=0.05)
	#fig.colorbar(iax, cax=cax, orientation='vertical')
	ax.set_xlabel("U")
	ax.set_ylabel(r"$\beta$")
	ax.set_title(r"$Z(\beta, U)$ hysteresis")
	plt.show()
	

