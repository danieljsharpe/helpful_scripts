'''
Python script to perform kMC for a given DTMC with associated reward values for edges, and calculate average total reward for transitions to absorbing node
'''

from __future__ import print_function
import numpy as np

'''
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
'''

'''
# renormalised network after eliminating node 1 of original chain
T = np.array([[0.810, 0.145, 0.045, 0.000],
              [0.180, 0.510, 0.160, 0.150],
              [0.020, 0.265, 0.715, 0.000],
              [0.000, 0.100, 0.000, 0.900]])

R = np.array([[ 0.0000000000,  0.178561284 , -0.8109302165,  0.0000000000],
              [-0.2557174   ,  0.0000000000,  0.05282132  , -0.4054651081],
              [ 0.8109302165, -0.785955730 ,  0.0000000000,  0.0000000000],
              [ 0.0000000000,  0.4054651081,  0.0000000000,  0.0000000000]])
'''

# renormalised network after eliminating node 2 of original chain
T = np.array([[0.64736842105, 0.2026315789, 0.15],
              [0.28026315789, 0.7197368421, 0.00],
              [0.10         , 0.00        , 0.90]])

R = np.array([[-0.01637209  , -0.18270332, -0.4054651081],
              [-0.68926471  ,  0.00000000,  0.0000000000],
              [ 0.4054651081,  0.00000000,  0.0000000000]])


seed=19 # seed for random number generator

b=2 # ID of initial node (indexed from 1)
a=3 # ID of target (absorbing) node

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
