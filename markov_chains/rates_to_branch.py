'''
Python script to read in network topology info from edge_conns.dat and edge_weights.dat files (DISCOTRESS format),
and construct the branching probability matrix and vector of mean waiting times
'''

from __future__ import print_function
import numpy as np

### INPUT
n=5 # no. of states in continuous-time Markov chain (CTMC)
e=6 # no. of edges (bidrectional transitions) in CTMC
do_committors = False # flag to also dump committor probability data from "committor_AB.dat" output file

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
K = P.copy() # transition rate matrix
for i in range(n):
    ksum = np.sum(P[i,:])
    tau[i] = 1./ksum # mean waiting time for transitions from i-th node
    K[i,i] = -1.*ksum
    P[i,:] *= tau[i]

print("\ntransition rate matrix:\n",K)
print("\nbranching probability matrix:\n",P)
print("\nmean waiting times vector:\n",tau)
print("\naction matrix (for shortest paths algorithm):\n",-1.*np.log(P))

# dump array info to files that can be read back in with pickle (e.g. via np.load())
K.dump("ratemtx.pkl")
P.dump("branchmtx.pkl")
tau.dump("meanwaittimes.pkl")

if not do_committors: quit()

q = np.zeros(n,dtype=float) # vector of committor probabilities
with open("committor_AB.dat","r") as q_f:
    i=0
    for line in q_f.readlines():
        q[i] = float(line.split()[0])
        i+=1
print("\ncommittor probabilities vector:\n",q)
q.dump("committors.pkl")
