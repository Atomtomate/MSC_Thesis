from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from distr import *
from matplotlib2tikz import save as tikz_save


plt.style.use('ggplot')
cmap = plt.get_cmap('jet_r')
class Sampler (object):
    def sample(self, N = 1):
        pass
        
    def expected_value(self):
        self.E_f = np.cumsum(self.f_sample*self.weights)/(np.cumsum(self.weights))
        
    def variance(self):
        n_step = np.arange(1,len(self.f_sample)+1)
        self.expected_value()
        print(self.f_sample*self.weights)
        #self.variance = np.cumsum( np.square(self.f_sample*self.weights - self.E_f) )/np.cumsum(self.weights) 
        self.variance = np.cumsum(self.f_sample*self.weights*self.f_sample*self.weights )/((np.cumsum(self.weights)-1.0))- np.square(self.E_f)
    
    def autocorr(self, maxt = 50):
        a = (self.f_sample*self.weights).copy()
        #m = a-np.mean(a)
        #norm = np.sum(m**2)
        #acorr = np.correlate(m, m, "same")/norm
        #self.acorr = acorr[int(len(acorr)/2):]
        #self.acorr = np.array([((a[:a.size-t] - self.E_f[-1])*(a[t:]- self.E_f[-1])).ravel().mean()/(self.variance[-1] * self.variance[-1]) for t in range(0,a.size)])
        #self.acorr = np.array([1]+[np.corrcoef(a[:-i], a[i:])[0,1] for i in range(1, int(a.size/2))])/(self.variance[-1])
        
    def run_sampler(self, N):
        res = self.sample(N)
        if(len(res.shape) > 1):
            self.f_sample = res[:,0]
            self.weights = res[:,1]
        else:
            self.f_sample = res
            self.weights = np.ones_like(self.f_sample)
        self.expected_value()
        self.variance()
        #self.autocorr()

    def show_results(self, true_distr = None, true_mean = None, true_variance = None, true_weight = None, title = "", save_tikz = False):
        fig, ax = plt.subplots()
        if true_distr is not None:
            x = np.linspace(np.min(self.f_sample*self.weights), np.max(self.f_sample*self.weights), len(self.f_sample)*55)
            ax.plot(x, true_distr(x), 'r-')
        ax.hist(self.f_sample, weights = self.weights, bins=np.arange(min(self.f_sample), max(self.f_sample) + 0.3, 0.3), normed = True, color=['royalblue'])
        ax.legend(['true distribution', 'sampled distribution'])
        if save_tikz:
              tikz_save(title.replace(" ", "") + "_hist.tex", figureheight = '8cm', figurewidth = '8cm')
        if true_weight is not None:
            fig2, ax_arr = plt.subplots(3, sharex=True)
        else:
            fig2, ax_arr = plt.subplots(2, sharex=True)
            
        #colors and stuff
        lw = 1.5
        
        #start plot
        ax_arr[0].set_xlim([-2,len(self.f_sample)+4])
        esl = ax_arr[0].plot(np.arange(1,len(self.f_sample)+1),np.abs(self.E_f), c='royalblue',lw=lw)#, c='orange'
        est = ax_arr[0].plot(np.arange(1,len(self.f_sample)+1), true_mean + 1.0/np.sqrt(np.arange(1,len(self.f_sample)+1)),ls='dotted', lw=lw)#, c='blue'
        if true_mean is not None:
            ax_arr[0].set_ylim([true_mean-max(1.0,true_mean)*0.8,true_mean+max(1.0,true_mean)*0.8])
        etl = ax_arr[0].axhline(true_mean, lw=(lw-0.5), c='red', ls = 'dashed')#
        ax_arr[0].legend(["$\mu_\mathrm{sampled}$", "$\mu_\mathrm{true}$", r'$\sqrt{\mathrm{samples}}^{-1}$'])
        ax_arr[0].set_xlabel("Samples")
        ax_arr[0].set_ylabel("$\mu$")
        
        nstep = np.arange(1,len(self.f_sample)+1)
        if true_variance is not None:
            ax_arr[1].set_ylim([true_variance-true_variance*0.8,true_variance+true_variance*0.8])
            vsl = ax_arr[1].plot(np.arange(1,len(self.variance)+1), self.variance, c='royalblue', lw=lw)
            vtl = ax_arr[1].axhline(true_variance, c='red', lw=(lw-0.5), ls = 'dashed')
            ax_arr[1].legend(["$\sigma^2_\mathrm{sampled}$", r'$\sigma^2_\mathrm{true}$'])
            ax_arr[1].set_xlabel("Samples")
            ax_arr[1].set_ylabel("Variance")
        
        if true_weight is not None:
            #ax_arr[2].set_ylim([true_weight-true_weight*0.2,true_weight+true_weight*0.2])
            wsl = ax_arr[3].plot(np.arange(1,len(self.weights)+1), np.cumsum(self.weights, dtype=float)/np.arange(1,len(self.weights)+1,dtype=float), c='royalblue', lw=lw)
            wtl = ax_arr[3].axhline(true_weight, c='red', lw=(lw-0.5), ls = 'dashed')
            ax_arr[3].legend(["$w_\mathrm{sampled}$", "$w_\mathrm{true}$"])
            ax_arr[3].set_xlabel("Samples")
            ax_arr[3].set_ylabel("Weight")
            
        #if save_tikz:
        #      tikz_save(title.replace(" ", "") + ".tex", figureheight = '8cm', figurewidth = '8cm')
    
        #fig3, ax3 = plt.subplots()
        #ax3.plot(np.arange(1,self.acorr.size + 1), self.acorr, c='royalblue', lw=lw)
          
        if save_tikz:
              tikz_save(title.replace(" ", "") + ".tex", figureheight = '8cm', figurewidth = '8cm')
        else:  
            plt.show()

        
class InversionSampler (Sampler):
    def __init__(self, right_inverse):
        self.right_inverse = right_inverse
        self.tries = 0
        self.samples = 0
        
    def sample(self, N = 1):
        self.tries += N
        self.samples += N
        return self.right_inverse(np.random.uniform(0.0, 1.0, N))
        
        
class RejectionSampler (Sampler):
    def __init__(self, p, envelope_distribution, scale = 1):
        self.p = p
        self.envelope_distribution = envelope_distribution
        self.scale = scale
        self.tries = 0
        self.samples = 0
        self.f_sample = []
        
    def sample(self, N = 1):
        for _ in range(N):
            while True:
                self.tries += 1
                x = self.envelope_distribution.rvs()   # draw a sample from envelope distr
                u = np.random.uniform(low=0, high=self.scale*self.envelope_distribution.pdf(x))
                if u < self.p(x):
                    self.samples += 1
                    self.f_sample.append(x)
                    break;
        return np.array(self.f_sample)
        

class ImportanceSampler (Sampler):
    def __init__(self, p, envelope_distribution):
        self.p = p
        self.envelope_distribution = envelope_distribution
        self.f_sample = []
    
    def sample(self, N = 1):
        for _ in range(N):
            x = self.envelope_distribution.rvs()
            #w = self.p(x)/self.envelope_distribution.pdf(x)
            w = np.exp( np.log(self.p(x))  - self.envelope_distribution.logpdf(x) )
            self.f_sample.append((x,w))
        return np.array(self.f_sample)


class MetropolisHastingsSampler (Sampler):
    def __init__(self, p, propose, acceptance, init_distr, burnin):
        self.p = p
        self.propose = propose
        self.acceptance = acceptance
        self.x = init_distr
        self.accepted = 0
        self.samples = 0
        self.f_sample = []
        self.burnin = burnin
        
    def sample(self, N = 1):
        for i in range(N):
            if i > self.burnin:
                self.samples += 1
            xp = self.propose(self.x)
            u = np.log(np.random.uniform())
            w = 1.0
            R = np.log(self.p(xp)) - np.log(self.p(self.x)) - np.log(self.acceptance(self.x,xp)) + np.log(self.acceptance(xp, self.x))
            #R = np.log(self.p(xp)) + np.log(self.acceptance(xp, self.x)) - np.log(self.p(self.x)) - np.log(self.acceptance(self.x,xp)) 
            
            if u < R:
                if i > self.burnin:
                    self.accepted += 1
                    self.f_sample.append((self.x,np.exp(w)))
                self.x = xp
        return np.array(self.f_sample)[self.burnin:,:]
    
        
def example_invSampl_normal_distribution():
    inv_sample_normal = InversionSampler(Phi_inv)
    inv_sample_normal.run_sampler(2000)
    inv_sample_normal.show_results(true_distr = phi, true_mean = 0.0, true_variance = 1.0, title = "Inversion Sampling" , save_tikz = True)
    
def example_rejSampl_normal_distribution():
    q = norm(loc = 1.2, scale = 2)
    scale = 5.05
    x = np.arange(-4,6,0.01)
    plt.plot(x, phi_test(x), 'b-', x, scale*q.pdf(x), 'g-')
    plt.legend(['Density p(x)', 'Envelope c*q(x)']);
    rej_sample_normal = RejectionSampler(phi_test, q, scale = scale)
    rej_sample_normal.run_sampler(100)
    rej_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3, title = "Rejection Sampling", save_tikz = False)
    print("tries: " + str(rej_sample_normal.tries) + ", samples: " + str(rej_sample_normal.samples) )
    
def example_impSampl_normal_distribution():
    q = norm(loc = 1.2, scale = 1.5)
    imp_sample_normal = ImportanceSampler(phi_test, q)
    imp_sample_normal.run_sampler(1000)
    imp_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3, true_weight = phi_test_norm(), title = "Importance Sampling", save_tikz = False)
    
def gaussian_propose(x):
    return norm(loc=x, scale=0.7).rvs()
def gaussian_step(x, xp):
    return norm(loc=x, scale=0.1).pdf(xp)
    
def example_MHSample_normal_distribution():
    MH_sample_normal = MetropolisHastingsSampler(phi_test, gaussian_propose, gaussian_step, 1.5, burnin = 1000)
    MH_sample_normal.run_sampler(20000)
    MH_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3, title = "Metropolis Sampling", save_tikz  = False)
    print("accepted samples: " + str(MH_sample_normal.accepted) + ", samples: " + str(MH_sample_normal.samples) + ", acceptance rate: " + str(MH_sample_normal.accepted/MH_sample_normal.samples))
    
example_invSampl_normal_distribution()
#example_rejSampl_normal_distribution()
#example_impSampl_normal_distribution()
#example_MHSample_normal_distribution()