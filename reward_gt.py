'''
Toy Python script to investigate the graph transformation (GT) algorithm to preserve the mean of edge-dependent quantities (rewards) by renormalisation

Daniel J. Sharpe
Nov 2020
'''

from __future__ import print_function
import numpy as np

''' compute the mean first passage rewards for the list of reward matrices Rmtxs using the iterative formulation of the GT algorithm '''
def calc_mfpr_iterative(n,T,Rmtxs):
    print("\nUsing iterative formulation of GT algorithm to compute MFPRs...")
    for x in range(n-2): # loop over nodes to be eliminated, leaving only the pair of nodes corresponding to the final two rows and columns in T
        print("\neliminating node:",x+1)
        for i in range(1,n-x):
            for j in range(1,n-x):
                Te = np.sum([T[0,k] for k in range(1,n-x)]) # indirectly evaluate factor of (1.-T[0,0])
                Tx = T[i,0]*T[0,j]/(1.-T[0,0]) # additional i->j transition probability in renormalisation
                # renormalise rewards
                if not (T[i,j]==0. and Tx==0.):
                    for Rmtx in Rmtxs: # loop over different reward matrices
                        dirp = T[i,j]/(T[i,j]+Tx) # conditional probability of direct transition
                        dr = dirp*Rmtx[i,j] # contribution from direct reward
                        ir = (1.-dirp)*(Rmtx[i,0]+Rmtx[0,j]+(((1./(1.-T[0,0]))-1.)*Rmtx[0,0])) # contribution from indirect reward
                        Rmtx[i,j] = dr+ir
                # renormalise transition probabilities
                T[i,j] += Tx # renormalise transition probabilities
        # eliminate state from transition and reward matrices
        T = T[1:,1:]
        Rmtxs = Rmtxs[:,1:,1:]
    rewards = np.zeros(np.shape(Rmtxs)[0],dtype=float)
    for k, Rmtx in enumerate(Rmtxs): # recall only two states remain, initial state (idx 0) and absorbing state (idx 1)
        Te = T[0,1]
        rewards[k] = (((1./Te)-1.)*Rmtx[0,0]) + Rmtx[0,1]
    return rewards

''' compute the mean first passage rewards for the list of reward matrices Rmtxs using the block formulation of the GT algorithm '''
def calc_mfpr_block(n,T,Rmtxs):
    print("\nUsing block formulation of GT algorithm to compute MFPRs...")
    rewards = np.zeros(np.shape(Rmtxs)[0],dtype=float)
    p0 = np.zeros(n-1,dtype=float) # initial occupation probability distribution vector
    p0[-1] = 1. # initial probability distribution localised at state corresponding to penultimate row and column of transition matrix
    G = np.linalg.inv(np.eye(n-1)-T[:-1,:-1]) # Green's (i.e. absorbing fundamental) matrix
    for k, Rmtx in enumerate(Rmtxs):
        Rvec = np.zeros(n-1,dtype=float) # average rewards for transitions from nodes
        for i in range(n-1): # recall final state is the absorbing state
            Rvec[i] = np.dot(T[i,:],Rmtx[i,:])
        rewards[k] = np.dot(np.dot(G,Rvec),p0)
    return rewards


if __name__=="__main__":

    ### INPUT
    # NB in the transition matrix, the final state is taken to be the absorbing state. The initial state can be chosen
    '''
    # specify transition probability matrix. The final state is the absorbing state and the script calculates the MFPT from a chosen initial state
    T = np.array([[0.50, 0.20, 0.15, 0.15, 0.00],
                  [0.15, 0.75, 0.10, 0.00, 0.00],
                  [0.20, 0.10, 0.45, 0.10, 0.15],
                  [0.05, 0.00, 0.25, 0.70, 0.00],
                  [0.00, 0.00, 0.10, 0.00, 0.90]])
    tau = 0.05 # lag time (or mean waiting time for linearised transition matrix of CTMC)
    tau_vec = [tau]*np.shape(T)[0]
    '''

    # alternatively, read in data from pickled branching probability matrix and vector of mean waiting times (NB entropy calc is then invalid)
    T = np.load("branchmtx.pkl")
    tau_vec = np.load("meanwaittimes.pkl")

    b = 1 # index of chosen initial state (from 1)
    iterative=True # compute mean first passage rewards using iterative formulation of GT algo (T/F) (if F, use block formulation)


    ### RUN
    n = np.shape(T)[0] # number of states
    # various reward matrices
    # time
    Rt = (np.repeat(tau_vec,n)).reshape((n,n))
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

    assert((b>0 and b<n))
    for i in range(n): assert abs(np.sum(T[i,:])-1.)<1.E-08
    # swap rows and columns of transition matrix so that initial state b corresponds to penultimate rows and cols
    idxs = [i for i in range(n-1) if i!=b-1]+[b-1,n-1]
    T = T[:,idxs]; T = T[idxs,:]
    for k, Rmtx in enumerate(Rmtxs):
        Rmtx = Rmtx[:,idxs]; Rmtx = Rmtx[idxs,:]
        Rmtxs[k] = Rmtx
    print("\ntransition matrix:\n",T)
    print("\nreward matrices (time,entropy,action):\n",Rmtxs[0],"\n\n",Rmtxs[1],"\n\n",Rmtxs[2])
    print("\ninitial state:",b,"final state:",n)

    if iterative: rewards = calc_mfpr_iterative(n,T,Rmtxs)
    else: rewards = calc_mfpr_block(n,T,Rmtxs)

print("\n\n")
print("Mean first passage time:\t",rewards[0])
print("Mean path entropy flow:\t\t",rewards[1])
print("Mean path action:\t\t",rewards[2])
