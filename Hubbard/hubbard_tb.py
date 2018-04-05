from __future__ import division
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as sp_linalg
import matplotlib.pyplot as plt
from scipy.sparse import dok_matrix, dia_matrix, csr_matrix, lil_matrix, coo_matrix
import seaborn as sns
from scipy.stats import norm

#sys.dont_write_bytecode = True

# === parameters ===
beta = 64.0
t    = 3.0
mu   = 0.0
dim  = 4
N    = [13,13,13,13]
Vol  = np.prod(N)
import time

Lp = N[1:]
    
def neib(pos, d):
    """ Compute neighbor of flattened index on a n-dimensional square lattice
    n - number of dimensions: n = L1*L2*...*LD*x0 + L2*L3*LD*x1 + ... + xD
    L - number of sites for each dimension [L0, L1, L2, ... LD]
    pos - flattened index of site
    d - [0,2n] direction of neighbor (>=n for negative direction)
    return - flattened index of neighbor
    """
    
    x  = [np.floor(pos/np.prod(Lp))]
    for i in range(1,dim):
        x.append(np.floor((pos - np.sum([np.prod(Lp[j:])*x[j] for j in range(i)]))/np.prod(Lp[i:])))
    dAbs = dim-d-1 if d >= 0 else d + dim 
    x[dAbs] = (x[dAbs] + (2*(d >= 0)-1))%N[dAbs]
    return np.sum([np.prod(Lp[i:])*np.array(x).astype(int)[i] for i in range(dim)])
    

def real_space():
    start_time = time.time()
    H = dok_matrix((Vol, Vol), dtype=np.float16)
    #H.setdiag(-mu)
    Hd =  [(i,neib(i, d)) for i in range(Vol) for d in range(-dim,dim)]
    for el in Hd:
        H[el[0],el[1]] = -t
    print("--- H generation took %s seconds ---" % (time.time() - start_time))
    H = H.tocsr() #.asfptype()
    start_time = time.time()
    vals = sp_linalg.eigsh(H,k=Vol-1,return_eigenvectors=False,tol=10e-1) 
    #vals, _ = linalg.eigh(H.todense())  
    print("--- Eigenproblem took %s seconds ---" % (time.time() - start_time))
    plt.plot(vals)
    plt.show()
    plt.hist(vals,bins=40)
    plt.show()
    sns.kdeplot(vals, bw=0.1)
    plt.show()
    np.save("vals_"+str(dim)+"D",vals)
    
    
    if dim == 1: 
        dat_1d = np.array([-mu - 2*t *np.cos(np.pi*i/N) for i in range(N-1)])
        plt.plot(dat_1d)
        plt.show()


#def k_space():
f, axarr = plt.subplots(2, 3, sharex= True)
for i,dim in enumerate([1,2,3,4,7,10]):
    a = np.ones(dim)
    ts = t/np.sqrt(2.0*dim)
    ax = axarr[0][i] if i < 3 else axarr[1][i-3]
    N = int(10e7)
    ek_sc = 2.0* ts * np.sum(np.cos(a[:,np.newaxis]*np.random.uniform(high=np.pi, size=(dim, N))), axis = 0)
    ek_bcc = 2.0* ts * np.prod(np.cos(0.5*a[:,np.newaxis]*np.random.uniform(high=np.pi, size=(dim, N))), axis = 0)
    dE = np.linspace(min(ek_sc),max(ek_sc),200)
    res = np.histogram(ek_sc, dE, density= True)[0]
    ax.plot(dE[1:], res, label=str(dim) + " dim.")
    if dim > 2:
        ax.plot(dE[1:], norm.pdf(dE[1:],0,t), label = r"$\mathcal{N}(0,1)$")
    ax.fill_between(dE[1:], 0, res, alpha = 0.5)
    ax.legend()
    ax.set_xlabel(r"$\epsilon$")
    ax.set_ylabel(r"$n(\epsilon)$")
    ax.set_yticks([])
plt.show()