from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


# parameters
beta = 64.0
U    = 2.5
mu   = 0.0

        
        
class HubbardAL:
    def __init__(self):
        self.s_up   = np.matrix([[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,1,0]])
        self.s_up_d = np.transpose(self.s_up)
        self.s_do   = np.matrix([[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,1,0,0]])
        self.s_do_d = np.transpose(self.s_do)
        self.n_up   = self.s_up_d * self.s_up
        self.n_do   = self.s_do_d * self.s_do
        
    
    def set_p(self, H, beta, U , mu):
        self.beta = beta
        self.U    = U
        self.mu   = mu
        self.H_diag = H
        self.Z      = np.trace(self.H_diag)

    
    def avg_n(self):
        avg_up = np.trace(self.H_diag * self.n_up)/self.Z
        avg_do = np.trace(self.H_diag * self.n_do)/self.Z
        var = np.trace(self.H_diag * ((self.n_up - self.n_do)**2) )/self.Z
        #avg_n  = np.trace(self.H_diag * (self.n_up+self.n_do))/self.Z#avg_up + avg_do
        avg_n = avg_up + avg_do
        return np.array([avg_n,var, avg_up, avg_do])


H = HubbardAL()
f, axarr = plt.subplots(2, 2)

for i,U in enumerate([4, 8]):
    x_list = np.linspace(-2.0*U,2.0*U, 1000)
    for beta in [0.5, 2, 32]:
        plt1 = []; plt2 = []
        for mu in x_list:
            H1 = np.matrix([[1.0,.0,.0,.0],[0.0,np.exp(beta*mu),0.0,0.0],[0.0,0.0,np.exp(beta*mu),0.0],[0.0,0.0,0.0,np.exp(-beta*(U-2.0*mu))]])
            H2 = np.matrix([[np.exp(-beta*U/4.0),.0,.0,.0],[.0,np.exp(beta*U/4.0 + beta*mu),.0,.0],[.0,.0,np.exp(beta*U/4.0 + beta*mu),.0],[.0,.0,.0,np.exp(-beta*U/4.0 + 2.0*beta*mu)]])
            H.set_p(H2, beta, U, mu)
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
        
    
mu = 0.0
f2, axarr2 = plt.subplots(2, 1)    
for U in [1, 2, 4, 8]:
    x_list = np.linspace(0.1,25, 2000)
    plt1 = []
    for beta in x_list:
        H2 = np.matrix([[np.exp(-beta*U/4.0),.0,.0,.0],[.0,np.exp(beta*U/4.0 + beta*mu),.0,.0],[.0,.0,np.exp(beta*U/4.0 + beta*mu),.0],[.0,.0,.0,np.exp(-beta*U/4.0 + 2.0*beta*mu)]])
        H.set_p(H2, beta, U, mu)
        res = H.avg_n()
        plt1.append(res[1])
    axarr2[0].semilogx(x_list, plt1, label = r"$U =$ " + str(U))
    axarr2[0].legend()
    axarr2[0].set_title(r"local moment at half-filling")
    axarr2[0].set_xlabel(r"$\beta$")
    axarr2[0].set_ylabel(r"$\Delta n $")
 
for beta in [0.5, 1, 3, 5]:
    x_list = np.linspace(0.1,100, 2000)
    plt1 = []
    for U in x_list:
        H2 = np.matrix([[np.exp(-beta*U/4.0),.0,.0,.0],[.0,np.exp(beta*U/4.0 + beta*mu),.0,.0],[.0,.0,np.exp(beta*U/4.0 + beta*mu),.0],[.0,.0,.0,np.exp(-beta*U/4.0 + 2.0*beta*mu)]])
        H.set_p(H2, beta, U, mu)
        res = H.avg_n()
        plt1.append(res[1])
    axarr2[1].semilogx(x_list, plt1, label = r"$\beta =$ " + str(beta))
    axarr2[1].legend()
    axarr2[1].set_title(r"local moment at half-filling")
    #axarr2[0].set_title(r"$U = $" + str(U))
    axarr2[1].set_xlabel(r"$U$")
    axarr2[1].set_ylabel(r"$\Delta n $")
                       
plt.show()
