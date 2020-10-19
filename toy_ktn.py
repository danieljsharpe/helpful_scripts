import numpy as np
from scipy.linalg import expm
from sys import argv

'''
A script to construct a toy kinetic transition network where the transition rates are given by an Arrhenius-type law (cf. harmonic
transition state theory). Such a Markov chain is guaranteed to satisfy the detailed balance condition (i.e. is reversible)
'''

write_output = int(argv[1]) # write output 1/0 (Y/N)

E = np.array([1.5,3.0,2.5,4.0,5.5]) # energies of discrete states
Eb = np.array([[0.0,4.0,5.0,5.0,0.0], # off-diagonal upper-triangular elements are energies of transition states (TSs)
               [0.0,0.0,4.5,0.0,0.0],
               [0.0,0.0,0.0,8.0,8.0],
               [0.0,0.0,0.0,0.0,0.0],
               [0.0,0.0,0.0,0.0,0.0]])

beta = 1. # effective inverse temperature
tau = 5. # lag time at which transition probability matrix is estimated

n = np.shape(E)[0]

K = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(i+1,n):
        if Eb[i,j]==0.: continue # zero entry in Eb indicates no connection between pair of discrete states
        assert np.array([Eb[i,j]-E[i]>0.,Eb[i,j]-E[j]>0.]).all() # TS energy must exceed energies of discrete states that it connects
        K[i,j] = np.exp((E[i]-Eb[i,j])*beta)
        K[j,i] = np.exp((E[j]-Eb[i,j])*beta)
for i in range(n): K[i,i] = -np.sum(K[i,:]) # ensures that *global* balance is satisfied (i.e. stationary distribution exists)
T = expm(K*tau)

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
print "\ntransition probability matrix:\n", T

if not write_output: quit()
evals_T, revecs_T = np.linalg.eig(T.T) # eigenvalues and right eigenvectors
idx = evals_T.argsort()[-1]
pi = revecs_T[:,idx]/np.sum(revecs_T[:,idx]) # stationary distribution
for i in range(n): assert abs(np.dot(pi,K)[i])<1.E-08
for i in range(n): assert abs(np.dot(pi,T)[i]-pi[i])<1.E-08


# an example problem to print
T=np.array([[ 0.5,   0.2,   0.15,  0.15,  0.  ],
            [ 0.15,  0.75,  0.1,   0.,    0.  ],
            [ 0.2,   0.1,   0.45,  0.1,   0.15],
            [ 0.05,  0.,    0.25,  0.7,   0.  ],
            [ 0.,    0.,    0.1,   0.,    0.9 ]])

K=np.array([[-10.,   4.,   3.,   3.,   0.],
            [  3.,  -5.,   2.,   0.,   0.],
            [  4.,   2., -11.,   2.,   3.],
            [  1.,   0.,   5.,  -6.,   0.],
            [  0.,   0.,   2.,   0.,  -2.]]) 

pi=np.array([0.15506773,0.20364316,0.19897244,0.14385801,0.29845866])
n = np.shape(T)[0]


# print DISCOTRESS input files for CTMC and DTMC
edge_weights_ctmc = open("edge_weights_ctmc.dat","w")
edge_conns_ctmc = open("edge_conns_ctmc.dat","w")
edge_weights_dtmc = open("edge_weights_dtmc.dat","w")
edge_conns_dtmc = open("edge_conns_dtmc.dat","w")
stat_prob_f = open("stat_prob.dat","w")
for i in range(n):
    stat_prob_f.write("%1.15f\n" % np.log(pi[i]))
    for j in range(i+1,n):
        if K[i,j]!=0.:
            edge_weights_ctmc.write("%1.15f  %1.15f\n" % (np.log(K[i,j]),np.log(K[j,i])))
            edge_conns_ctmc.write("%4i %4i\n" % (i+1,j+1))
        if T[i,j]!=0.:
            edge_weights_dtmc.write("%1.15f  %1.15f\n" % (T[i,j],T[j,i]))
            edge_conns_dtmc.write("%4i %4i\n" % (i+1,j+1))
edge_weights_ctmc.close()
edge_conns_ctmc.close()
edge_weights_dtmc.close()
edge_conns_dtmc.close()
stat_prob_f.close()
