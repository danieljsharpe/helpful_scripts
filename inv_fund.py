''' toy python script to demonstrate various computations (e.g. inversion to obtain fundamental matrix for
    both irreducible and absorbing Markov chains), starting from a discrete-time Markov chain and also showing
    the continuous-time case'''

import numpy as np
from math import sqrt

# discrete-time stochastic transition matrix - the last state should be absorbing (no outgoing transiitons)
P = np.array([[0.5,0.2,0.15,0.15,0.],[0.15,0.75,0.1,0.,0.],[0.2,0.1,0.45,0.1,0.15],
              [0.05,0.,0.25,0.7,0.],[0.,0.,0.1,0.,0.9]],dtype=float)

'''
# Kemeny and Snell's textbook example
p = 2./3.
q = 1.-p
M = np.array([[1.,-p,0.],[-q,1.,-p],[0.,-q,1.]],dtype=float)
Q = M.copy()
'''

for i in range(np.shape(P)[0]): # check stochasticity
    assert abs(np.sum(P[i,:])-1.) < 1.E-12

Q = P[:-1,:-1] # substochastic matrix of nonabsorbing states (i.e. Q-matrix of transition matrix in canonical form)
tau = 0.05 # lag time of DTMC

### INVESTIGATE FUNDAMENTAL MATRIX OF ABSORBING DTMC

''' in the following, note that the element n_{ij} of the fundamental matrix N is equal to the expected number of times
    node j is visited, given that the process starts in node i, prior to absorption. This does not include the fact that
    the process starts in node i, so the row sums are equal to the number of steps when starting in that state '''

M = np.eye(np.shape(Q)[0],dtype=float)-Q # Markovian kernel

N = np.linalg.inv(M) # fundamental matrix of absorbing Markov chain
Nvar = np.dot(N,2.*np.diag(np.diagonal(N))-np.eye(np.shape(Q)[0]))-(N*N)
H = np.dot(N-np.eye(np.shape(Q)[0]),np.diag(np.array([1./x for x in np.diagonal(N)])))

l = np.dot(N,np.ones(np.shape(Q)[0])) # vector of expected path lengths
lvar =  np.dot((2.*N)-np.eye(np.shape(Q)[0]),l)-(l*l)  # vector of variances associated with MFPTs

print "\ndiscrete-time transition probability matrix:\n", P
print "\nfundamental matrix of absorbing chain (mean numbers of node visits):\n", N
print "\nvariances in numbers of node visits:\n", Nvar
print "\nvisitation probability matrix:\n", H
print "\n"
for i in range(np.shape(Q)[0]):
    print "expected time to absorption (MFPT) / variance thereof, when starting in state %i:    %.6f   /   %.6f" \
        % (i+1,l[i]*tau,lvar[i]*tau)

### EIGENDECOMPOSITION OF DISCRETE-TIME CHAIN

evals, revecs = np.linalg.eig(P.T) # calculate right eigenvectors of the irreducible stochastic matrix
#print "\nright eigenvectors before sorting:\n", revecs
revecs = np.array([revecs[i,:] for i in evals.argsort()[::-1]])
evals, levecs = np.linalg.eig(P) # calculate left eigenvectors of the irreducible stochastic matrix
#print "\nleft eigenvectors before sorting:\n", levecs
levecs = np.array([levecs[:,i] for i in evals.argsort()[::-1]]).T
print "\neigenvalues\n", evals[evals.argsort()[::-1]]
print "\nright eigenvectors:\n", revecs[:,:]
#print "\nleft eigenvectors:\n", levecs[:,:]
assert abs(evals[evals.argsort()[-1]]-1.)<1.E-14 # there should be a single eigenvalue equal to unity
pi = revecs[:,0]/np.sum(revecs[:,0]) # equilibrium occupation probabilities
for i in range(np.shape(P)[0]):
    for j in range(i,np.shape(P)[0]):
#        print "i:", i+1, "j:", j+1, "\t", np.dot(revecs[:,i],revecs[:,j])
        print "i:", i+1, "j:", j+1, "\t", np.dot(revecs[:,i],levecs[:,j])
# compute normalisation factors
rnorm = np.zeros(np.shape(P)[0],dtype=float) # normalisation factors associated with right eigenvectors

tmp_arr_right = np.zeros((np.shape(P)[0],np.shape(P)[0]),dtype=float) # diagonal elements required for normalisation of right eigenvectors
tmp_arr_left = np.zeros((np.shape(P)[0],np.shape(P)[0]),dtype=float) # diagonal elements required for normalisation of left eigenvectors
for i in range(np.shape(P)[0]):
    for j in range(np.shape(P)[0]):
        for k in range(np.shape(P)[0]):
            tmp_arr_right[i,j] += revecs[k,i]*revecs[k,j]/pi[k]
            tmp_arr_left[i,j] += levecs[k,i]*levecs[k,j]*pi[k]
# normalise
for i in range(np.shape(P)[0]):
    revecs[:,i] *= 1./sqrt(tmp_arr_right[i,i])
    levecs[:,i] *= 1./sqrt(tmp_arr_left[i,i])
# '''
# check orthonormality of both left and right eigenvectors
print "\ntmp_arr_right:\n", tmp_arr_right
quit()
print "\ntmp_arr_left:\n", tmp_arr_left
for i in range(np.shape(P)[0]):
    print abs(np.sum(revecs[:,i])-(lambda i: 1. if i==0 else 0.)(i))
#    assert abs(np.sum(revecs[:,i])-(lambda i: 1. if i==0 else 0.)(i))<1.E-10
#    assert abs(np.dot(pi,levecs[:,i])-(lambda i: 1. if i==0 else 0.)(i))<1.E-10
    for j in range(i,np.shape(P)[0]):
        if i==j: continue
#        assert abs(tmp_arr_right[i,j])<1.E-10
#        assert abs(tmp_arr_left[i,j])<1.E-10
# '''

MFPT = np.zeros((np.shape(P)[0],np.shape(P)[0]),dtype=float)
for i in range(np.shape(P)[0]):
    for j in range(np.shape(P)[0]):
        for k in range(1,np.shape(P)[0]):
            MFPT[i,j] += (1./(pi[i]*np.log(evals[k])))*revecs[k,i]*(levecs[k,i]-levecs[k,j])

print "\nnormalised right eigenvectors:\n", revecs
print "\nnormalised left eigenvectors:\n", levecs
print "\nstationary distribution:\n", pi
quit()
print "\nmatrix of MFPTs:\n", MFPT

### INVESTIGATE FUNDAMENTAL MATRIX OF IRREDUCIBLE DTMC

A = 1. # note that this is actually the group inverse, Kemeny and Snell's fundamental matrix Z
