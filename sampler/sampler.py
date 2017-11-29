from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from distr import *


class Sampler (object):
    def sample(self, N = 1):
        pass
        
    def expected_value(self):
        self.E_f = np.cumsum(self.f_sample*self.weights)/(np.cumsum(self.weights))
        
    def variance(self):
        n_step = np.arange(1,len(self.f_sample)+1)
        self.expected_value()
        #self.variance = np.cumsum( np.square(self.f_sample*self.weights - self.E_f) )/np.cumsum(self.weights) 
        self.variance = np.cumsum(self.f_sample*self.f_sample*self.weights )/(np.cumsum(self.weights)) - np.square(self.E_f)
        
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

    def show_results(self, true_distr = None, true_mean = None, true_variance = None, true_weight = None, title = ""):
        fig, ax = plt.subplots()
        if true_distr is not None:
            x = np.linspace(np.min(self.f_sample*self.weights), np.max(self.f_sample*self.weights), len(self.f_sample)*10)
            ax.plot(x, true_distr(x), 'r-')
        ax.hist(self.f_sample, weights = self.weights, normed = True)
        ax.legend(['true distribution', 'sampled distribution'])
        
        if true_weight is not None:
            fig2, ax_arr = plt.subplots(3, sharex=True)
        else:
            fig2, ax_arr = plt.subplots(2, sharex=True)
        
        esl = ax_arr[0].plot(np.arange(1,len(self.f_sample)+1), self.E_f)
        if true_mean is not None:
            ax_arr[0].set_ylim([true_mean-true_mean*0.4,true_mean+true_mean*0.4])
        etl = ax_arr[0].axhline(true_mean, color = 'r')
        ax_arr[0].legend(['E_sampled', 'E_true'])
        
        nstep = np.arange(1,len(self.f_sample)+1)
        if true_variance is not None:
            ax_arr[1].set_ylim([true_variance-true_variance*0.4,true_variance+true_variance*0.4])
            vsl = ax_arr[1].plot(np.arange(1,len(self.variance)+1), self.variance)
            vtl = ax_arr[1].axhline(true_variance, color = 'r', ls = 'dashed')
            ax_arr[1].legend(['Var_sampled', 'Var_true'])
        
        if true_weight is not None:
            #ax_arr[2].set_ylim([true_weight-true_weight*0.2,true_weight+true_weight*0.2])
            wsl = ax_arr[2].plot(np.arange(1,len(self.weights)+1), np.cumsum(self.weights, dtype=float)/np.arange(1,len(self.weights)+1,dtype=float))
            wtl = ax_arr[2].axhline(true_weight, color = 'r', ls = 'dotted')
            ax_arr[2].legend(['weight_sampled', 'weight_true'])
            
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
    inv_sample_normal.run_sampler(1000)
    inv_sample_normal.show_results(true_distr = phi, true_mean = 0.0, true_variance = 1.0)
    
def example_rejSampl_normal_distribution():
    q = norm(loc = 1.2, scale = 2)
    scale = 5.05
    x = np.arange(-4,6,0.01)
    plt.plot(x, phi_test(x), 'b-', x, scale*q.pdf(x), 'g-')
    plt.legend(['Density p(x)', 'Envelope c*q(x)']);
    rej_sample_normal = RejectionSampler(phi_test, q, scale = scale)
    rej_sample_normal.run_sampler(6000)
    rej_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3)
    print("tries: " + str(rej_sample_normal.tries) + ", samples: " + str(rej_sample_normal.samples) )
    
def example_impSampl_normal_distribution():
    q = norm(loc = 1.2, scale = 1.5)
    imp_sample_normal = ImportanceSampler(phi_test, q)
    imp_sample_normal.run_sampler(3000)
    imp_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3, true_weight = phi_test_norm())
    
def gaussian_propose(x):
    return norm(loc=x, scale=1.8).rvs()
def gaussian_step(x, xp):
    return norm(loc=x, scale=0.4).pdf(xp)
    
def example_MHSample_normal_distribution():
    MH_sample_normal = MetropolisHastingsSampler(phi_test, gaussian_propose, gaussian_step, 3.0, burnin = 1000)
    MH_sample_normal.run_sampler(8000)
    MH_sample_normal.show_results(true_distr = phi_test_true, true_mean = 1.2, true_variance = 1.3*1.3)
    print("accepted samples: " + str(MH_sample_normal.accepted) + ", samples: " + str(MH_sample_normal.samples) + ", acceptance rate: " + str(MH_sample_normal.accepted/MH_sample_normal.samples))
    
example_invSampl_normal_distribution()
example_rejSampl_normal_distribution()
example_impSampl_normal_distribution()
example_MHSample_normal_distribution()