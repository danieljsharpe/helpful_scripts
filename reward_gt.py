'''
Toy Python script to investigate the graph transformation (GT) algorithm to preserve the mean of edge-dependent quantities (rewards) by renormalisation

Daniel J. Sharpe
Nov 2020
'''

from __future__ import print_function
import numpy as np

# transition probability matrix. The final state is the absorbing state and the script calculates the MFPT from the initial state
T = np.array([[0.50, 0.20, 0.15, 0.15, 0.00],
              [0.15, 0.75, 0.10, 0.00, 0.00],
              [0.20, 0.10, 0.45, 0.10, 0.15],
              [0.05, 0.00, 0.25, 0.70, 0.00],
              [0.00, 0.00, 0.10, 0.00, 0.90]])

n = np.shape(T)[0] # number of states

# various reward matrices
# time
tau = 0.05
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
# gather reward matrices in list
Rp = Rp.reshape((n,n))
Rmtcs = np.array((Rt,Rs,Rp))

print("original T:\n",T)
idxd = [0,3,2,1,4]
T  = T[:,idxd]; T = T[idxd,:]
print("new T:\n",T)
print("\n\nidxd:\n",idxd)
#quit()

### RUN ###
print("\ntransition matrix:\n",T)
print("\nreward matrices (time,entropy,action):\n",Rmtcs[0],"\n\n",Rmtcs[1],"\n\n",Rmtcs[2])
# swap the final state to correspond to the first row and first column of the transition matrix
idx = [i for i in range(-1,n-1)]
T = T[:,idx]; T = T[idx,:]
for Rmtx in Rmtcs:
    Rmtx = Rmtx[:,idx]; Rmtx = Rmtx[idx,:]
print("\nT:\n",T)
# graph transformation
for x in range(n-1,1,-1): # loop over nodes to be eliminated except absorbing and initial nodes
    print("\n\neliminating node:",x)
    # renormalise rewards
    for i in range(x):
        for j in range(x):
            print("  ",idxd[i]+1,idxd[j]+1)
            Tx = T[i,x]*T[x,j]/(1.-T[x,x]) # additional i->j transition probability in renormalisation
            if not (T[i,j]==0. and Tx==0.):
                for k, Rmtx in enumerate(Rmtcs): # loop over different reward matrices
                    if k!=1: continue
#                    print("    ",T[i,j],"    ",Tx)
#                    print(T[i,j]/(T[i,j]+Tx))
                    dirp = T[i,j]/(T[i,j]+Tx) # conditional probability of direct transition
                    dr = dirp*Rmtx[i,j] # contribution from direct reward
                    ir = (1.-dirp)*(Rmtx[i,x]+Rmtx[x,j]+(((1./(1.-T[x,x]))-1.)*Rmtx[x,x])) # contribution from indirect reward
                    print("      dir mtx elem:\t",T[i,j])
                    print("      cprob direct:\t",dirp)
                    print("      direct reward:\t",Rmtx[i,j])
                    print("      indirect reward:\t",Rmtx[i,x]+Rmtx[x,j]+(((1./(1.-T[x,x]))-1.)*Rmtx[x,x]))
                    print("      avg reward:\t",dr+ir)
#                    print("      contribution from direct reward:\t",dr)
#                    print("      contribution from indirect reward:\t",ir)
                    Rmtx[i,j] = dr+ir
    # renormalise transition probabilities
    for i in range(x):
        for j in range(x):
            Tx = T[i,x]*T[x,j]/(1.-T[x,x])
            T[i,j] += Tx # renormalise transition probabilities
    # eliminate state from transition and reward matrices
    T = T[:-1,:-1]
    for Rmtx in Rmtcs: Rmtx = Rmtx[:-1,:-1]
    print("\nT:\n",T,"\nR:\n",Rmtcs[1])
    quit()

rewards = np.zeros(n,dtype=float)
for i, Rmtx in enumerate(Rmtcs): # recall only two states remain, initial state (idx 1) and absorbing state (idx 0)
    rewards[i] = ((1./(1.-T[1,1])-1.)*Rmtx[1,1]) + T[1,0]
print("Mean first passage time:\t",rewards[0])
print("Mean path entropy flow:\t\t",rewards[1])
print("Mean path action:\t\t",rewards[2])
