'''
Python script to read in network topology info from edge_conns.dat and edge_weights.dat files (DISCOTRESS format),
and construct the branching probability matrix and vector of mean waiting times
'''

from __future__ import print_function
import numpy as np

### INPUT
n=8 # no. of states in continuous-time Markov chain (CTMC)
e=16 # no. of edges (bidrectional transitions) in CTMC

### RUN
conns = np.zeros((e,2),dtype=int) # list of connections in CTMC
with open("edge_conns.dat","r") as conns_f:
    i=0
    for line in conns_f.readlines():
        conns[i,0] = int(line.split()[0])-1
        conns[i,1] = int(line.split()[1])-1
        i+=1
wts = np.zeros((e,2),dtype=float)
with open("edge_weights.dat","r") as wts_f: # note that file contains ln transition rates
    i=0
    for line in wts_f.readlines():
        wts[i,0] = np.exp(float(line.split()[0]))
        wts[i,1] = np.exp(float(line.split()[1]))
        i+=1

P = np.zeros((n,n),dtype=float) # branching probability matrix
tau = np.zeros(n,dtype=float) # vector of mean waiting times
for conn, wt in zip(conns,wts):
    P[conn[0],conn[1]] = wt[0]
    P[conn[1],conn[0]] = wt[1]
for i in range(n):
    tau[i] = 1./np.sum(P[i,:])
    P[i,:] *= tau[i]

print("\nbranching probability matrix:\n",P)
print("\nmean waiting times vector:\n",tau)

# dump array info to files that can be read back in with pickle (e.g. via np.load())
P.dump("branchmtx.pkl")
tau.dump("meanwaittimes.pkl")
