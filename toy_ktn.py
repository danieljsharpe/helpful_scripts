import numpy as np
from scipy.linalg import expm

'''
A script to construct a toy kinetic transition network where the transition rates are given by an Arrhenius-type law (cf. harmonic
transition state theory). Such a Markov chain is guaranteed to satisfy the detailed balance condition (i.e. is reversible)
'''

E = np.array([1.5,3.0,2.5,4.0,5.5]) # energies of discrete states
Eb = np.array([[0.0,4.0,5.0,5.0,0.0], # off-diagonal upper-triangular elements are energies of transition states (TSs)
               [0.0,0.0,4.5,0.0,0.0],
               [0.0,0.0,0.0,8.0,8.0],
               [0.0,0.0,0.0,0.0,0.0],
               [0.0,0.0,0.0,0.0,0.0]])

T = 1. # effective temperature
tau = 5. # lag time at which transition probability matrix is estimated

n = np.shape(E)[0]

K = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(i+1,n):
        if Eb[i,j]==0.: continue # zero entry in Eb indicates no connection between pair of discrete states
        assert np.array([Eb[i,j]-E[i]>0.,Eb[i,j]-E[j]>0.]).all() # TS energy must exceed energies of discrete states that it connects
        K[i,j] = np.exp((E[i]-Eb[i,j])/T)
        K[j,i] = np.exp((E[j]-Eb[i,j])/T)
for i in range(n): K[i,i] = -np.sum(K[i,:]) # ensures that *global* balance is satisfied (i.e. stationary distribution exists)

# check that *detailed* balance condition holds
evals, revecs = np.linalg.eig(K.T)
idx = evals.argsort()[-1]
assert abs(evals[idx])<1.E-08 # stationary process is associated with zero eigenvalue
pi = revecs[:,idx]/np.sum(revecs[:,idx])
for i in range(n):
    for j in range(i+1,n):
        assert abs((K[i,j]*pi[i])-(K[j,i]*pi[j]))<1.E-08

np.set_printoptions(precision=12)
print "\ntransition rate matrix:\n", K
print "\ntransition probability matrix:\n", expm(K*tau)
