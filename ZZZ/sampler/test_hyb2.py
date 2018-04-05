from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

U	= 4.
beta	= 30.
tsteps	= 120
mu 	= U/2.
D	= 1.0
t       = D/2.
avg_n   = 5.2
g0mS	= 1.j*np.loadtxt("0_G0_Guess_MF.out", skiprows=1, usecols=(0,2))
g0mS_r	= np.loadtxt("0_G0_Guess_MF.out", skiprows=1, usecols=(0,2))

def SE_tail(U, ns, wn):
	return U*(ns - 0.5) + U*U*ns*(1.-ns)/(1.j*wn)

def mf(n, b):
	return (2.0*n+1)*np.pi/b

def tail(wn, cm1, c1, c3):
	res = cm1*wn - c1/wn - c3/(wn*wn*wn)
	return res

def tailFit(inp, orders, usable):
	xdata = np.array([mf(u, beta) for u in range(usable,len(inp))])
	ydata = inp[usable:]
	popt, pcov = curve_fit(tail, xdata, ydata)
	return popt, pcov

def trafo(inp, tau, mmm = [0.,0.,0.]):
	res = 0. + 0.j
	for el in inp:
		res += ((el[1]) * np.exp(-el[0]*tau)- mmm[1]/(el[0]))   #+  1.j*moments[2]/(el[0]*el[0]*el[0]) # - 1.0/(1.0j*el[0])
	res = (res - 0.5*mmm[1])/beta# +  0.*moments[2]*tau*(beta-tau)*0.25
	return res

def symToFull(inp):
	res = np.array([(-el[0],-el[1]) for el in inp[::-1]])[:-1]
	res = np.vstack((res, inp))
	return res

def GtoD(inp, c1 = 0.25):
	res = []
	for el in inp:
		#input file has low resolution (only first moment fit). 
		#dirty fix, possible because output is generated from Bethe lattice:
		if abs(el[1] - 1.0/(el[1])) < 1E-5:
			res.append((el[0], c1/el[0]))
		else:
			res.append((el[0], el[0] - 1.0/(el[1]))) 
	return np.array(res)

def fitTest(inp, fitHyb = False):
	if not np.isreal(inp).any():
		inp = np.imag(inp)
	u_start = 50
	u = int(inp.shape[0]/2)+u_start
	x1 = int(inp.shape[0]) - 100
	
	par, cov = tailFit(inp[:,1], 3, u)
	print(par)
	print(cov)
	x = mf(np.arange(u_start,u_start+u + 0.8*len(inp[:,1])), beta)
	y = tail(x, *par)
	f, ax = plt.subplots()
	ax.plot(inp[:,0],inp[:,1], label=r"original data", marker='.')
	ax.plot(x,y,label = r"fitted tail")
	if fitHyb:
		par2 = [0., D*D/4., 0.]
		y2 = tail(x, *par2)
		ax.plot(x,y2,label = r"tail from DOS moments")
	ax.set_xlabel(r"$i \omega_n$")
	ax.set_ylabel(r"$G(i \omega_n)$")
	ax.legend()
	#sub_axes = plt.axes([.6, .3, .25, .25])
	#xD = mf(np.arange(x1,x1+60), beta)
	#yD = tail(xD, *par)
	#sub_axes.plot(inp[x1:x1+60,0], inp[x1:x1+60,1], marker='.')
	#sub_axes.plot(xD, inp[x1:x1+60,1]-yD)
	#sub_axes.set_ylim(-0.0011738535,-0.0011738525)
	#sub_axes.set_ylim(min(inp[x1,1],y[x1])-(inp[x1,1]-y[x1]), max(inp[x1,1],y[x1])+(inp[x1,1]-y[x1]))
	#sub_axes.set_xlim(inp[x1+1,0],inp[x1+1,0]+1)
	#sub_axes.set_yticks([0,yD[0]])
	#sub_axes.set_yticklabels([str(inp[x1,1]),str(yD[0])])
	#sub_axes.set_ylabel(r"$3\cdot 10^{-7} - 3.5\cdot 10^{-7}$")
	plt.show()



def plot_hyb():
	#hybmS  = np.loadtxt("0_Hyb_MF.out", skiprows=1, usecols=(0,2))
	#hybtS  = np.array([trafo(hybmS, t*beta/tsteps) for t in range(tsteps)])
	g0mS 	= 1.j*np.loadtxt("0_G0_Guess_MF.out", skiprows=1, usecols=(0,2))
	g0m 	= symToFull(g0mS)
	DmS     = GtoD(g0mS)
	Dm      = GtoD(g0m)
	g0tS  = np.array([trafo(g0mS, t*beta/(tsteps)) for t in range(1,tsteps)])
	g0t  = np.array([trafo(g0m, t*beta/(tsteps)) for t in range(1,tsteps)])
	par, cov = tailFit(np.imag(g0mS[:,1]), 0, 10)
	print par
	print cov
	g0tTail = np.array([trafo(g0m, t*beta/(tsteps), [0., 1.,0.]) for t in range(1,tsteps)])
	Dt   = np.array([trafo(Dm, t*beta/(tsteps)) for t in range(1,tsteps)])
	par2, cov2 = tailFit(np.imag(DmS[:,1]), 0, 10)
	print par2
	print cov2
	DtTailE  = np.array([trafo(DmS, t*beta/(tsteps),[0., D*D/4.0,0.]) for t in range(1,tsteps)])
	DtTail   = np.array([trafo(DmS, t*beta/(tsteps), par2) for t in range(1,tsteps)])
	f, axarr = plt.subplots(1, sharex=True)
	axarr.set_title(r"$G(i \omega_N)$")
	axarr.plot(np.imag(g0mS[:,0]),np.imag(g0mS[:,1]), marker='.')
	f, axarr = plt.subplots(1, sharex=True)
	axarr.set_title(r"$\Delta (i \omega_N)$")
	axarr.plot(np.imag(DmS[:,0]),np.imag(DmS[:,1]), 'bo', markersize=1.1)
	axarr.set_xlabel(r"$\omega_n$")
	plt.show()
	#f2, axarr2 = plt.subplots(1, sharex=True)
	plt.plot(np.linspace(1,beta,tsteps-1),np.real(g0t), marker='.', label=r"$G_0 (\tau)$")
	plt.plot(np.linspace(1,beta,tsteps-1),np.real(g0tTail), marker='.', label=r"Tail corrected $G_0(\tau)$")
	plt.legend()
	#axarr2.plot(np.linspace(1,beta,tsteps-1),np.real(Dt), marker='.', label=r"$\Delta (\tau)$")
	#axarr2[1].plot(np.linspace(1,beta,tsteps-1),np.real(DtTail), marker='.', label=r"Tail corrected $\Delta (\tau)$")
	#axarr2[1].legend()
	plt.xticks([1,beta], ['0',r"$\beta$"])
	plt.xlabel(r"$\tau$")
	plt.show()

#plot_hyb()
g0m 	= symToFull(g0mS)
DmS     = GtoD(g0mS)
Dm      = GtoD(g0m)

#wnl = np.array([mf(i, beta) for i in range(-30,31)])
#seTail = SE_tail(U, avg_n, wnl)
#plt.plot(wnl, np.imag(se_tail)); plt.show()
