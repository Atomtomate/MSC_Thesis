from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np

# Normal distribution, CDF and right inverse
def phi(x):
    return norm.pdf(x)

def Phi(x):
    return norm.cdf(x)

def Phi_inv(x):
    return norm.ppf(x)
    
def plot_norm():
    x = np.arange(-5, 5, 0.01)
    plt.figure()
    plt.plot(x, phi(x), 'g-')
    plt.plot(x, Phi(x), 'b-')
    
    u = np.arange(0, 1, 0.01)
    plt.figure()
    plt.plot(u, Phi_inv(u), 'b-')
 
mean = 1.2
std = 1.3
def phi_test(x, mean = 1.2, std = 1.3 ):
    return np.exp((-(x-mean)**2.0)/(2.0*std*std) )
    
def phi_test_norm():
    return np.sqrt(2.0*np.pi*std*std)
    
def phi_test_true(x):
    return np.exp((-(x-mean)**2.0)/(2.0*std*std) )/phi_test_norm()