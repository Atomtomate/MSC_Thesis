from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpmath import *
import matplotlib.mlab as mlab
#import scipy.misc as smisc

# parameters
beta = 64.0
U    = 2.5
mu   = 0.0
eta  = 0.025

# normal np.matrix([[1.0,.0,.0,.0],[0.0,np.exp(beta*mu),0.0,0.0],[0.0,0.0,np.exp(beta*mu),0.0],[0.0,0.0,0.0,np.exp(-beta*(U-2.0*mu))]])
# phs np.matrix([[np.exp(-beta*U/4.0),.0,.0,.0],[.0,np.exp(beta*U/4.0 + beta*mu),.0,.0],[.0,.0,np.exp(beta*U/4.0 + beta*mu),.0],[.0,.0,.0,np.exp(-beta*U/4.0 + 2.0*beta*mu)]]) 
        
class HubbardAL:
    def __init__(self, beta, U , mu, phs = True):
        self.s_up   = np.array([[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,1,0]])
        self.s_up_d = np.transpose(self.s_up)
        self.s_do   = np.array([[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,1,0,0]])
        self.s_do_d = np.transpose(self.s_do)
        self.n_up   = np.matmul(self.s_up_d, self.s_up)
        self.n_do   = np.matmul(self.s_do_d, self.s_do)
        self.beta = beta
        self.U    = U
        self.mu   = mu
        self.phs  = phs
        
    
    def set_p(self, beta, U , mu):
        self.beta = beta
        self.U    = U
        self.mu   = mu
        self.H_diag = self.get_H(beta)
        self.Z      = np.trace(self.H_diag,axis1=1, axis2=2) if len(self.H_diag.shape)>2 else np.trace(self.H_diag)
        
    def get_H(self, x):
        z = np.zeros_like(x)
        one = np.ones_like(x)
        if self.phs:
            res = np.array([np.array([np.exp(-x*self.U/4.0),z,z,z]),\
                            np.array([z,np.exp(x*(self.U/4.0 + self.mu)),z,z]),\
                            np.array([z,z,np.exp(x*(self.U/4.0 + self.mu)),z]),\
                            np.array([z,z,z,np.exp(-x*(self.U/4.0 - 2.0*self.mu))])])
        else:
            res = np.array([np.array([one,z,z,z]),\
                            np.array([z,np.exp(x*self.mu),z,z]),\
                            np.array([z,z,np.exp(x*self.mu),z]),\
                            np.array([z,z,z,np.exp(-x*(self.U-2.0*self.mu))])])
        return np.rollaxis(res,2,0) if len(res.shape)>2 else res
        
    def gf(self, tau):
        H1 = self.get_H(self.beta - np.array(tau))
        H2 = self.get_H(tau) 
        gt = -np.trace(np.matmul(np.matmul(np.matmul(H1,self.s_up),H2),self.s_up_d),axis1=1, axis2=2)/self.Z
        #gt = (np.exp(self.mu*tau) + np.exp(self.beta * self.mu - tau*(self.U - self.mu)))/self.Z
        gw = np.array([(1+ np.exp(self.beta*self.mu))/(self.mu + complex(0.,1.)*(2.*n + 1)*np.pi/self.beta) +\
                      (np.exp(-self.beta*(self.U - 2.*self.mu)) + np.exp(self.beta*self.mu))/(self.mu -self.U + complex(0.,1.)*(2.*n + 1)*np.pi/self.beta) for n in range(-20,20)]/self.Z)
        _, _, avg_u, avg_d = self.avg_n()
        se_u = self.U*avg_d + selfU*self.U*(avg_d*(1.-avg_d))/(omega + self.mu - self.u*(1.-avg_d))
        #gw_d = np.array([(1.0 - avg_u)/(complex(0.,1.)*(2.*n + 1)*np.pi/self.beta + self.mu) + avg_u/(complex(0.,1.)*(2.*n + 1)*np.pi/self.beta + self.mu - self.U)  for n in range(-20,20) ])
        return np.array([gt, gw])
        #print(H1.shape)
     
    def selfEnergy(self, w):
        _, _, avg_u, avg_d = self.avg_n()
        se_u = self.U*avg_d + selfU*self.U*(avg_d*(1.-avg_d))/(w + self.mu - self.u*(1.-avg_d))
        se_d = self.U*avg_u + selfU*self.U*(avg_u*(1.-avg_u))/(w + self.mu - self.u*(1.-avg_u))   
        return [se_u, se_d]
    
    def avg_n(self):
        avg_up = np.trace(self.H_diag * self.n_up)/self.Z
        avg_do = np.trace(self.H_diag * self.n_do)/self.Z
        var = np.trace(self.H_diag * ((self.n_up - self.n_do)**2) )/self.Z
        #avg_n  = np.trace(self.H_diag * (self.n_up+self.n_do))/self.Z#avg_up + avg_do
        avg_n = avg_up + avg_do
        return np.array([avg_n,var, avg_up, avg_do])


H = HubbardAL(beta, U , mu, phs = True)

def plot_rho():
    f, axarr = plt.subplots(2, 2)
    
    for i,U in enumerate([4, 8]):
        x_list = np.linspace(-2.0*U,2.0*U, 1000)
        for beta in [0.5, 2, 32]:
            plt1 = []; plt2 = []
            for mu in x_list:
                H.set_p(beta, U, mu)
                res = H.avg_n()
                plt1.append(res[0])
                plt2.append(res[1])
            axarr[0][i].plot(x_list, plt1, label = r"$\beta =$ " + str(beta))
            axarr[0][i].legend()
            axarr[0][i].set_title(r"$U = $" + str(U))
            if i == 0:
                axarr[0][i].set_ylabel(r"$\langle n \rangle$")
                axarr[1][i].set_ylabel(r"$\Delta n $")
            axarr[1][i].plot(x_list, plt2, label = r"$\beta =$ " + str(beta))
            #axarr[i][0].legend()
            #axarr[i][0].set_title(r"$U = $" + str(U))
            axarr[1][i].set_xlabel(r"$\mu$")
            axarr[0][i].get_shared_x_axes().join(axarr[0][i], axarr[1][i])
        
def plot_local_moment():  
    mu = 0.0
    f, axarr = plt.subplots(1, 2, sharey = True)    
    for U in [1, 2, 4, 8]:
        x_list = np.linspace(0.1,25, 2000)
        plt1 = []
        for beta in x_list:
            H.set_p(beta, U, mu)
            res = H.avg_n()
            plt1.append(res[1])
        axarr[0].semilogx(x_list, plt1, label = r"$U =$ " + str(U))
        axarr[0].legend()
        #axarr2[0].set_title(r"local moment at half-filling")
        axarr[0].set_xlabel(r"$\beta$")
        axarr[0].set_ylabel(r"$\Delta n $")
    
    for beta in [0.4, 1, 8]:
        x_list = np.linspace(0.1,100, 2000)
        plt1 = []
        for U in x_list:
            H.set_p(beta, U, mu)
            res = H.avg_n()
            plt1.append(res[1])
        axarr[1].semilogx(x_list, plt1, label = r"$\beta =$ " + str(beta))
        axarr[1].legend()
        #axarr2[1].set_title(r"local moment at half-filling")
        #axarr2[0].set_title(r"$U = $" + str(U))
        axarr[1].set_xlabel(r"$U$")
        #axarr2[1].set_ylabel(r"$\Delta n $")               
        plt.show()
        
        
def plot_gf():
    x = np.linspace(0,1,400)
    f, axarr = plt.subplots(1,2)
    for beta in [1,8, 32]:
        tau = np.linspace(0,beta,1000)
        H.set_p(beta, U , 0.0)
        g_tau, gw = H.gf(tau)
        axarr[0].plot(tau/beta, g_tau ,label =  r"$\beta =$ " + str(beta))
        axarr[1].plot(range(-int(len(gw)/2),int(len(gw)/2)) ,np.imag(gw),  marker = '.', linestyle= 'none', label =  r"$\beta =$ " + str(beta))
    axarr[0].legend()
    axarr[0].set_xlabel(r"$\tau/\beta$")
    axarr[0].set_ylabel(r"$G(\tau)$", labelpad=-20)
    axarr[0].locator_params(nticks=2)
    axarr[0].set_xticks([0.0,0.5, 1.0])
    axarr[0].set_yticks([-0.5, 0.])
    axarr[1].legend()
    axarr[1].set_xlabel(r"$\omega_n$")
    axarr[1].set_ylabel(r"$Im[G(i \omega)]$")
    axarr[1].set_yticks([])
    plt.show()
    
    
    x = np.linspace(-3.*U/4.,3.*U/4.,1000)
    y = mlab.normpdf(x, -U/2., 0.025)+mlab.normpdf(x, U/2., eta)
    plt.plot(x,y)
    plt.fill(x, y, alpha=0.5)
    plt.ylabel(r"$A(\omega)$")
    plt.xlabel(r"$\omega$", labelpad=-10)
    plt.xticks([-U/2, U/2],['-U/2', 'U/2'])
    plt.yticks([])
    plt.show()
    
#plot_gf()

x = np.linspace(-U,U,1000)
plt.plot(x, np.real( (U*U/(4*(x+complex(0,1.)*eta)))), label = r"$Re[ \Sigma (\omega) ]$")
plt.plot(x, np.imag( (U*U/(4*(x+complex(0,1.)*eta)))), label = r"$Im[ \Sigma (\omega) ]$")
plt.xticks([0],['0'])
plt.yticks([0],['0'])
plt.ylabel(r"$\Sigma(\omega)$")
plt.xlabel(r"$\omega$", labelpad=-10)
plt.legend()
plt.show()