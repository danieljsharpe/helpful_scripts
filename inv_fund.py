''' toy python script to invert fundamental matrix '''

import numpy as np


# Q-matrix of transition matrix in canonical form
Q = np.array([[0.5,0.2,0.15,0.15,0.],[0.15,0.75,0.1,0.,0.],[0.2,0.1,0.45,0.1,0.15],
              [0.05,0.,0.25,0.7,0.],[0.,0.,0.1,0.,0.9]],dtype=float)

M = np.eye(np.shape(Q)[0],dtype=float)-Q # Markovian kernel

N = np.linalg.inv(M) # fundamental matrix of absorbing Markov chain

#print "\nfundamental matrix:\n", N

''' in the following, note that the element n_{ij} of the fundamental matrix N is equal to the expected number of times
    node j is visited, given that the process starts in node i, prior to absorption. This does not include the fact that
    the process starts in node i, so the row sums are equal to the number of steps when starting in that state '''
#print "\n"
#for i in range(np.shape(N)[0]):
#    print "number of steps before absorption, when starting in state %i:    %.6f" % (i,np.sum(N[:,i]))


evals, evecs = np.linalg.eig(Q.T) # calculate right eigenvectors of Q
print "\neigenvalues:\n", evals
print "\neigenvectors:\n", evecs

print "\nstat distribn:\n", evecs[:,0]/np.sum(evecs[:,0])
