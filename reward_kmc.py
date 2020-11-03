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

b=4 # ID of initial node (indexed from 1)
a=5 # ID of target (absorbing) node (indexed from 1)

seed=19 # seed for random number generator

npaths = 100000 # no. of paths to simulate

tau=0.05 # lag time

# various reward matrices
n = np.shape(T)[0] # number of nodes
# time
Rt = (np.repeat(tau,n*n)).reshape((n,n))
# entropy flow
Rs = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(i,n):
        if i==j or T[i,j]==0.: continue # assuming all edges are bidirectional
        Rs[i,j] = -np.log(T[i,j]/T[j,i])
        Rs[j,i] = -np.log(T[j,i]/T[i,j])
# path action
Rp = -1.*np.log(T.flatten())
nan_idx = np.argwhere(np.isinf(Rp))
Rp[nan_idx] = 0.
# gather reward matrices in array
Rp = Rp.reshape((n,n))
Rmtxs = np.array([Rt,Rs,Rp])
m = np.shape(Rmtxs)[0] # number of reward matrices
for k in range(m): print(Rmtxs[k])

### run simulation

for i in range(n): assert abs(np.sum(T[i,:])-1.)<1.E-08
assert (b>=1 and b<=n)
assert (a>=1 and a<=n)
np.random.seed(seed)
Rvals = np.zeros((m,npaths),dtype=float) # list of reward values for paths (time,entropy,action)
# main simulation loop
for i in range(npaths):
    x = b-1 # x denotes current node
    Rpath = np.zeros(m,dtype=float) # rewards for current iteration
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
        for k in range(m): # rewards for transition
            Rpath[k] += Rmtxs[k,x,y]
        x = y
    for k in range(m):
        Rvals[k,i] = Rpath[k]

Ravg = np.zeros(m,dtype=float) # averages of rewards for first passage path ensemble
for k in range(m): Ravg[k] = np.sum(Rvals[k])/float(npaths) # average reward
print("\nMean first passage time (MFPT):\t",Ravg[0])
print("\nAverage path entropy flow:\t",Ravg[1])
print("\nAverage path action:\t\t",Ravg[2])
