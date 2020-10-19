from __future__ import print_function
import numpy as np
from scipy import linalg
from math import sqrt

#'''
# note: this rate matrix satisfies global, but not detailed, balance
K = np.array([[-5.,2.0,2.0,1.0],
              [1.5,-6.,1.5,3.0],
              [3.0,4.0,-9.,2.0],
              [0.5,0.5,1.0,-2.]])
tau = 5.E-03 # lag time
#'''

'''
# note: this rate matrix satisfies detailed balance
K = np.array([[-0.142479765469,  0.082084998624,    0.030197383422,  0.030197383422,  0.            ],
              [ 0.367879441171,  -0.59100960132,    0.223130160148,  0.,              0.            ],
              [ 0.082084998624,  0.135335283237,   -0.225593824737,  0.004086771438,  0.004086771438],
              [ 0.367879441171,  0.            ,    0.018315638889, -0.38619508006,   0.            ],
              [ 0.,              0.,                0.082084998624,  0.,             -0.082084998624]])
tau = 5. # lag time
'''

n = np.shape(K)[0] # number of states

evals_K, revecs_K = np.linalg.eig(K.T)
idx = evals_K.argsort()[-1]
pi = revecs_K[:,idx]/np.sum(revecs_K[:,idx])
evals_K = evals_K[evals_K.argsort()[::-1]]
D = np.diag([1./pi[i] for i in range(n)])

assert abs(np.sum(pi)-1.)<1.E-08
for i in range(n): assert abs(np.dot(pi,K)[i])<1.E-08 # check global balance condition
print("\ntransition rate matrix:\n",K)
print("\nstationary distribution:\n",pi)
print("\neigenvalues:\n",evals_K)
# check if detailed balance condition satisfied (not necessary for fundamental matrices, Kemeny constant etc. to be valid)
reversible=True
for i in range(n):
    for j in range(n):
        if abs((K[i,j]*pi[i])-(K[j,i]*pi[j]))>1.E-08: reversible=False
print("\ndetailed balance condition satisifed?",reversible)

Z = np.linalg.inv(np.outer(np.ones(n),pi)-K)-np.outer(np.ones(n),pi) # fundamental matrix (some freedom in choice, here gives Kemeny constant directly as trace)
MFPT = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        MFPT[i,j] = (1./pi[j])*(Z[j,j]-Z[i,j]) # note MFPT matrix for CTMC has diagonal elems equal to zero
for i in range(n): assert abs(np.dot(pi,MFPT[i,:])-np.trace(Z))<1.E-08 # check alt definition of Kemeny constant
assert abs((np.sum([1./abs(x) for x in evals_K[1:]]))-np.trace(Z))<1.E-08 # Kemeny constant from CTMC eigenvalues
# check np.dot(MFPT,pi) is a left eigenvector of the transition rate matrix associated with eigenvalue equal to zero
kp = np.dot(K,np.dot(MFPT,pi))
for i in range(n): assert abs(kp[i])<1.E-08
# MFPT matrix from CTMC should recover precisely the same rate matrix as the original
KK = np.dot(D-np.ones((n,n)),np.linalg.inv(MFPT))
for kk_elem, k_elem in zip(KK.flatten(),K.flatten()): assert abs(kk_elem-k_elem)<1.E-08

print("\nfundamental matrix from CTMC:\n",Z)
print("\nKemeny constant from CTMC:\n", np.trace(Z))
print("\nMFPT matrix from CTMC:\n",MFPT)

### DISCRETE-TIME MARKOV CHAIN

T = linalg.expm(K*tau) # discrete-time transition probability matrix
evals_T, revecs_T = np.linalg.eig(T.T)
idx = evals_T.argsort()[-1]
piT = revecs_T[:,idx]/np.sum(revecs_T[:,idx])
evals_T = evals_T[evals_T.argsort()[::-1]]
for i in range(n): assert abs(np.dot(pi,T)[i]-pi[i])<1.E-08 # check global balance condition
for i in range(n): assert abs(piT[i]-pi[i])<1.E-08 # check stationary distribution is identical to CTMC result
print("\ntransition probability matrix:\n",T)

ZT = np.linalg.inv(np.eye(n)-T+np.outer(np.ones(n),pi)) # Kemeny and Snell's fundamental matrix
MFPT_T = np.dot(np.eye(n)-ZT+np.dot(np.ones((n,n)),np.diag(np.diagonal(ZT))),D) # note MFPT matrix for DTMC has non-zero diagonal elems
for i in range(n): assert abs(np.dot(pi,MFPT_T[i,:])-np.trace(ZT))<1.E-08 # check alt definition of Kemeny constant (in terms of MFPT matrix and stat distribn)
#assert abs((1.+np.sum([1./(1.-x) for x in evals_T[1:]]))-np.trace(ZT))<1.E-08 # Kemeny constant from DTMC eigenvalues
# check np.dot(MFPT_T,pi) is a left eigenvector of the transition probability matrix associated with eigenvalue equal to unity
tp1 = np.dot(MFPT_T,pi)
tp2 = np.dot(T,np.dot(MFPT_T,pi))
for i in range(1,n): assert abs(tp1[0]-tp1[i])<1.E-08 # all elems of tp1 (and tp2) should be equal
for i in range(n): assert abs(tp1[i]-tp2[i])<1.E-08

print("\nKemeny constant from DTMC:\n",np.trace(ZT)*tau)
print("\nMFPT matrix from DTMC:\n", MFPT_T*tau)

# recover transition rate matrix from the MFPT matrix computed from the DTMC
KT = np.dot(D-np.ones((n,n)),np.linalg.inv((MFPT_T-np.diag(np.diagonal(MFPT_T)))*tau))
for i in range(n): assert abs(np.sum(KT[i,:]))<1.E-08
''' note that there is a many-to-one mapping of CTMCs to DTMCs. Check that the recovered
    CTMC corresponds to the DTMC that also corresponds to the original CTMC. Further note: the
    various CTMCs may also correspond to different branching probability matrices '''
TT = linalg.expm(KT*tau)
for tt_elem, t_elem in zip(TT.flatten(),T.flatten()): print("\nWARNING: recovered transition matrix does not match original"); break # assert abs(tt_elem-t_elem)<1.E-08
print("\nrecovered transition matrix:\n",TT)
print("\nrecovered transition rate matrix:\n",KT)

# recover transition probability matrix from the MFPT matrix computed from the CTMC
pi_mfpt = np.dot((MFPT+np.diag(D))/tau,np.diag(pi))
pi_arr = np.outer(np.ones(n),pi)
TK = np.dot(pi_arr-pi_mfpt,np.linalg.inv(np.eye(n)-pi_mfpt))
print("\ntransition probability matrix from CTMC MFPTs:\n",TK)
for i in range(n): assert abs(np.sum(TK[i,:])-1.)<1.E-08 # check row-stochasticity
#for tk_elem, t_elem in zip(TK.flatten(),T.flatten()): assert abs(tk_elem-t_elem)<1.E-08

# BRANCHING PROBABILITY MATRICES

# note that the branching probability matrices for the different transition rate matrices are not necessarily the same
BK = K.copy()-np.diag(np.diagonal(K)) # branching probability matrix from original transition rate matrix
BT = KT.copy()-np.diag(np.diagonal(KT)) # branching probability matrix from recovered transition rate matrix
for i in range(n):
    BK[i,:] *= 1./-K[i,i]
    BT[i,:] *= 1./-KT[i,i]
print("\nBranching matrix from original rate matrix:\n",BK)
print("\nBranching matrix from recovered rate matrix:\n",BT)

# COMPUTE MFPTs FROM THE ABSORBING FUNDAMENTAL MATRIX

''' function to compute NxN matrix of MFPTs from transition probability matrix by considering a set of N (N-1)x(N-1) substochastic
    reducible transition matrices. tau_vec is a vector of waiting times; for a DTMC, each entry is equal to the lag time, for a CTMC,
    the j-th entry is the mean waiting time for the j-th node '''
def mfpt_mtx_reduc(Pm,tau_vec):
    n = np.shape(Pm)[0]
    idcs = np.array([i for i in range(n)]) # keep track of indices of nodes when swapping rows and columns of transition matrix
    MFPT_mtx = np.zeros((n,n),dtype=float)
    for i in range(n): # cycle through target states
        Pm_copy = Pm.copy()
        Pm_copy[[-1,i],:] = Pm_copy[[i,-1],:] # swap i-th row with final row, the latter is considered absorbing
        Pm_copy[:,[-1,i]] = Pm_copy[:,[i,-1]] # swap i-th col with final col, the latter is considered absorbing
        tau_vec_copy = tau_vec.copy()
        tau_vec_copy[[-1,i]] = tau_vec[[i,-1]] # swap elems of vector of mean waiting times
        idcs_copy = idcs.copy()
        idcs_copy[[-1,i]] = idcs[[i,-1]]
        N = np.linalg.inv(np.eye(n-1)-Pm_copy[:-1,:-1])
        MFPT_vals = np.dot(N,tau_vec_copy[:-1])
        d=0
        for j in range(n): # write the computed MFPTs to the array to be returned
            if j==i:
                d=1; continue
            MFPT_mtx[idcs_copy[j-d],i] = MFPT_vals[j-d]
    return MFPT_mtx
            
MFPT_BK = mfpt_mtx_reduc(BK,np.array([-1./x for x in np.diagonal(K)]))
for mfpt_bk, mfpt_ctmc in zip(MFPT_BK.flatten(),MFPT.flatten()): assert abs(mfpt_bk-mfpt_ctmc)<1.E-04
MFPT_BT = mfpt_mtx_reduc(BT,np.array([-1./x for x in np.diagonal(KT)]))
MFPT_BT += D*tau
for mfpt_bt, mfpt_dtmc in zip(MFPT_BT.flatten(),MFPT_T.flatten()*tau): assert abs(mfpt_bt-mfpt_dtmc)<1.E-04

ZBK = np.linalg.inv(np.eye(n)-BK+np.outer(np.ones(n),pi)) # fundamental matrix from branching probability matrix

MFPT_Tr = mfpt_mtx_reduc(T,np.array([tau for i in range(n)]))
MFPT_Tr += D*tau
for mfpt_tr, mfpt_dtmc in zip(MFPT_Tr.flatten(),MFPT_T.flatten()*tau): assert abs(mfpt_tr-mfpt_dtmc)<1.E-04


# QUANTITIES PERTAINING TO DETAILED BALANCE CONDITION

dS = np.zeros((n,n),dtype=float) # matrix of entropy flow along edges (as multiple of Boltzmann constant)
for i in range(n):
    for j in range(i+1,n):
        if K[i,j]==0. and K[j,i]==0.: continue # nodes are not connected in either direction
        dS[i,j] = -np.log(K[i,j]/K[j,i])
        dS[j,i] = -np.log(K[j,i]/K[i,j])
print("\nentropy flow matrix:\n",dS)

if not reversible: quit()
''' define the symmetrized rate matrix. Note that the symmetrized form does not have rows that sum to zero;
    nonetheless, the symmetrized rate matrix has the same eigenvalues as the original rate matrix '''
K_sym = np.zeros((n,n),dtype=float)+np.diag(np.diagonal(K))
for i in range(n):
    for j in range(n):
        if i==j: continue
        K_sym[i,j] = sqrt(pi[i]/pi[j])*K[i,j]
evals_sym, evecs_sym = np.linalg.eig(K_sym.T)
evecs_sym = np.array([evecs_sym[:,i] for i in evals_sym.argsort()[::-1]])
evals_sym = evals_sym[evals_sym.argsort()[::-1]]
for i in range(n): assert abs(evals_sym[i]-evals_K[i])<1.E-08 # symmetrised and original rate matrices have same eigenvalues
for i in range(n): assert abs((pi[i]/sqrt(pi[i]))-evecs_sym[0,i])<1.E-08 # test relation for right eigenvector of reversible stochastic matrix on first eigvec
for i in range(n): assert abs((1.*sqrt(pi[i]))-evecs_sym[0,i])<1.E-08 # test relation for left eigevector of reversible stochastic matrix on first eigvec
for i in range(n): # eigenvectors of symmetrised rate matrix are orthonormal
    for j in range(n):
        assert abs(np.dot(evecs_sym[i,:],evecs_sym[j,:])-(lambda i,j: 1. if i==j else 0.)(i,j))<1.E-08
print("\neigenvectors of symmetrised rate matrix:\n",evecs_sym)
