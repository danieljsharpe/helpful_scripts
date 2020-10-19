'''
Python script to perform kMC for a given DTMC with associated reward values for edges, and calculate average total reward for transitions to absorbing node
'''

from __future__ import print_function
import numpy as np

# transition probability matrix
T = np.array([[0.50, 0.20, 0.15, 0.15, 0.00],
              [0.15, 0.75, 0.10, 0.00, 0.00],
              [0.20, 0.10, 0.45, 0.10, 0.15],
              [0.05, 0.00, 0.25, 0.70, 0.00],
              [0.00, 0.00, 0.10, 0.00, 0.90]])

# reward matrix
R =  np.array([[ 0.0000000000, -0.2876820725,  0.2876820725, -1.098612289,  0.0000000000],
               [ 0.2876820725,  0.0000000000,  0.0000000000,  0.000000000,  0.0000000000],
               [-0.2876820725,  0.0000000000,  0.0000000000,  0.916290732, -0.4054651081],
               [ 1.0986122890,  0.0000000000, -0.9162907319,  0.000000000,  0.0000000000],
               [ 0.0000000000,  0.0000000000,  0.4054651081,  0.000000000,  0.0000000000]])

seed=19 # seed for random number generator

b=4 # ID of initial node (indexed from 1)
a=5 # ID of target (absorbing) node

npaths = 100000 # no. of paths to simulate

tau=0.05 # lag time

### run simulation

n = np.shape(T)[0] # number of nodes
for i in range(n): assert abs(np.sum(T[i,:])-1.)<1.E-08
assert n==np.shape(R)[0]
assert (b>=1 and b<=n)
assert (a>=1 and a<=n)
np.random.seed(seed)
rvals = np.zeros(npaths,dtype=float) # list of reward values for paths
tvals = np.zeros(npaths,dtype=int) # list of MFPTs for paths
# main simulation loop
for i in range(npaths):
    x = b-1 # x denotes current node
    r = 0. # reward of current iteration
    n = 0 # number of steps in path
#    print("path no:",i)
    while x!=a-1:
#        print("  x:",x+1)
        s = np.random.rand()
#        print("    s:",s)
        tsum = 0.
        for z, t in enumerate(T[x,:]):
            tsum += t
            if s < tsum:
                y = z # y denotes next node
                break
#        print ("    y:",y+1)
        r += R[x,y] # reward for transition
        n += 1
        x = y
    rvals[i] = r
    tvals[i] = n

ravg = np.sum(rvals)/float(npaths) # average reward
print("\naverage reward:\t",ravg)
mfpt = np.sum(tvals.astype(float)*tau)/float(npaths) # mean first passage time (MFPT)
print("\nMFPT:\n",mfpt)
