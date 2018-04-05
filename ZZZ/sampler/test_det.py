from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

np.set_printoptions(precision = 4, linewidth = 300)

def Mdec(M, k = 1):
	Rp = M[-k:,:-k]
	Qp = M[:-k,-k:]
	Sp = M[-k:,-k:]
	Pp = M[:-k,:-k]
	Mp = Pp - Qp*Rp/(Sp)
	return Mp

def Minc(M, R, Q, S):
	t = M.dot(Q); tt = R.dot(M)
	if(S.shape[0] > 1):
		Sp = inv(S - R.dot(t))
		Rp = -Sp.dot(tt)
		Pp = M + t.dot((Sp.dot(tt)))
	else:
		Sp = np.array([1./(S - R.dot(t))])
		Rp = -Sp*tt
		Pp = M + t.dot(Sp*tt)
	Qp = -t.dot(Sp)
	tmpT = np.hstack((Pp, Qp))
	tmpR = np.hstack((Rp, Sp))
	Mp   = np.vstack((tmpT, tmpR))
	return Mp

M4 = np.array([[1, 2, 3, 3.5],[ 4, 5, 6, 6.5],[ 7, 8, 9, 9.5],[10,11,12,12.5]])

R4 = np.array([1.1, 2.1, 3.1, 4.1])
Q4 = np.array([[1.2], [2.2], [3.2], [4.2]])
S4 = np.array([1.3])

R5 = np.array([1.4, 2.4, 3.4, 4.4, 5.4])
Q5 = np.array([[1.5], [2.5], [3.5], [4.5], [5.5]])
S5 = np.array([1.6])

R46 = np.array([[1.1, 2.1, 3.1, 4.1], [1.4, 2.4, 3.4, 4.4]])
Q46 = np.array([[1.2, 1.5],[2.2, 2.5], [3.2, 3.5], [4.2, 4.5]])
S46 = np.array([[1.3, 5.5],[5.4, 1.6]])

M5  = Minc(M4, R4, Q4, S4)
M6  = Minc(M5, R5, Q5, S5)
M6p = Minc(M4, R46, Q46, S46)
M4p = Mdec(M5)
np.testing.assert_almost_equal(M4, M4p)
np.testing.assert_almost_equal(M6, M6p)
M4pp = Mdec(M6, 2)
np.testing.assert_almost_equal(M4, M4pp)
