import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn')
sns.set_style("whitegrid")

test = np.array([70.        ,  2.05      , -0.09191591, -0.24332052, -0.35910384, -0.44235926, -0.49888241, -0.53474862])
wnlist = np.array([1j*(2.*n+1.)*np.pi/test[0] for n in range(len(test)-2)])
test2 = np.array([70.        ,  3.05      , -0.74226073, -1.84369238, -2.3795746 , -2.53100616, -2.489231  , -2.36816348])
wnlist2 = np.array([1j*(2.*n+1.)*np.pi/test2[0] for n in range(len(test)-2)])

res = np.polyfit(wnlist, 1j*test[2:], 4).imag
res2 = np.polyfit(wnlist2, 1j*test2[2:], 4).imag
res = np.poly1d(res)
res2= np.poly1d(res2)
resd = np.polyder(res)
res2d = np.polyder(res2)
x = np.linspace(0.,0.1,100)
plt.figure()
plt.title(r"$\beta = $" + str(test[0]))
plt.plot(x, res(x), label="U=2.05")
plt.plot(x, res2(x),label="U=3.05")
plt.xlabel(r"$\omega$")
plt.ylabel(r"$Im[\Sigma(\omega)]$")
plt.legend()
plt.figure()
plt.title(r"$\beta = $" + str(test[0]))
plt.plot(x, resd(x), label="U=2.05")
plt.plot(x, res2d(x),label="U=3.05")
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\frac{\partial\, Im[\Sigma(\omega)]}{\partial \omega}$")
plt.legend()
plt.show()
