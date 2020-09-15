import numpy as np
from scipy import linalg

K = np.array([[-5.,2.0,2.0,1.0],
              [1.5,-6.,1.5,3.0],
              [3.0,4.0,-9.,2.0],
              [0.5,0.5,1.0,-2.]])
n = np.shape(K)[0]

evals_K, revecs_K = np.linalg.eig(K.T)
pi = revecs_K[:,0]/np.sum(revecs_K[:,0])
D = np.diag([1./pi[i] for i in range(n)])

for i in range(n): assert abs(np.dot(pi,K)[i])<1.E-08
print "\ntransition rate matrix:\n", K
print "\nstationary distribution:\n", pi

Z = np.linalg.inv(np.outer(np.ones(n),pi)-K)-np.outer(np.ones(n),pi) # fundamental matrix (some freedom in choice, here gives Kemeny constant directly as trace)
MFPT = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        MFPT[i,j] = (1./pi[j])*(Z[j,j]-Z[i,j])
for i in range(n): assert abs(np.dot(pi,MFPT[i,:])-np.trace(Z))<1.E-08 # check alt definition of Kemeny constant
assert abs((np.sum([1./abs(x) for x in evals_K[1:]]))-np.trace(Z))<1.E-08 # Kemeny constant from CTMC eigenvalues
# check np.dot(MFPT,pi) is a left eigenvector of the transition rate matrix associated with eigenvalue equal to zero
kp = np.dot(K,np.dot(MFPT,pi))
for i in range(n): assert abs(kp[i])<1.E-08
KK = np.dot(D-np.ones((n,n)),np.linalg.inv(MFPT)) # this should recover precisely the same rate matrix as the original
for kk_elem, k_elem in zip(KK.flatten(),K.flatten()): assert abs(kk_elem-k_elem)<1.E-08

print "\nfundamental matrix from CTMC:\n", Z
print "\nKemeny constant from CTMC:\n", np.trace(Z)
print "\nMFPT matrix from CTMC:\n", MFPT

### DISCRETE-TIME MARKOV CHAIN

tau = 0.05 # lag time
T = linalg.expm(K*tau) # discrete-time transition probability matrix
print "\ntransition probability matrix:\n", T
evals_T, revecs_T = np.linalg.eig(T.T)
piT = revecs_T[:,0]/np.sum(revecs_T[:,0])
for i in range(n): assert abs(piT[i]-pi[i])<1.E-08 # check stationary distribution is identical to CTMC result

ZT = np.linalg.inv(np.eye(n)-T+np.outer(np.ones(n),pi)) # Kemeny and Snell's fundamental matrix
MFPT_T = np.dot(np.eye(n)-ZT+np.dot(np.ones((n,n)),np.diag(np.diagonal(ZT))),D)
for i in range(n): assert abs(np.dot(pi,MFPT_T[i,:])-np.trace(ZT))<1.E-08 # check alt definition of Kemeny constant (in terms of MFPT matrix and stat distribn)
assert abs((1.+np.sum([1./(1.-x) for x in evals_T[1:]]))-np.trace(ZT))<1.E-08 # Kemeny constant from DTMC eigenvalues
MFPT_nodiag = MFPT_T-np.diag(np.diagonal(MFPT_T))
# check np.dot(MFPT_T,pi) is a left eigenvector of the transition probability matrix associated with eigenvalue equal to unity
tp1 = np.dot(MFPT_T,pi)
tp2 = np.dot(T,np.dot(MFPT_T,pi))
for i in range(1,n): assert abs(tp1[0]-tp1[i])<1.E-08 # all elems of tp1 (and tp2) should be equal
for i in range(n): assert abs(tp1[i]-tp2[i])<1.E-08

print "\nKemeny constant from DTMC:\n", np.trace(ZT)*tau
print "\nMFPT matrix from DTMC:\n", MFPT_T*tau

KT = np.dot(D-np.ones((n,n)),np.linalg.inv((MFPT_T-np.diag(np.diagonal(MFPT_T)))*tau))
for i in range(n): assert abs(np.sum(KT[i,:]))<1.E-08
''' note that there is a many-to-one mapping of CTMCs to DTMCs. Check that the recovered
    CTMC corresponds to the DTMC that also corresponds to the original CTMC. Further note: the
    various CTMCs may also correspond to different branching probability matrices '''
TT = linalg.expm(K*tau)
for tt_elem, t_elem in zip(TT.flatten(),T.flatten()): assert abs(tt_elem-t_elem)<1.E-08
print "\nrecovered transition rate matrix:\n", KT

# note that the branching probability matrices for the different transition rate matrices are not necessarily the same
BK = K.copy()-np.diag(np.diagonal(K)) # branching probability matrix from original transition rate matrix
BT = KT.copy()-np.diag(np.diagonal(KT)) # branching probability matrix from recovered transition rate matrix
for i in range(n):
    BK[i,:] *= 1./-K[i,i]
    BT[i,:] *= 1./-KT[i,i]
print "\nBranching matrix from original rate matrix:\n", BK
print "\nBranching matrix from recovered rate matrix:\n", BT

# recover transition probability matrix from MFPTs
T1 = np.dot(np.outer(np.ones(n),pi)-np.dot(MFPT_T,np.diag(pi)),np.linalg.inv(np.eye(n)-np.dot(MFPT_T,np.diag(pi))))
for i in range(n): assert abs(np.sum(T1[i,:])-1.)<1.E-08 # check row-stochasticity
for t1_elem, t_elem in zip(T1.flatten(),T.flatten()): assert abs(t1_elem-t_elem)<1.E-08
print "\nT1:\n", T1
pi_arr = np.outer(np.ones(n),pi)
DT_arr = np.dot((MFPT/tau)+D,np.diag(pi)) # could instead use MFPT_T to get back expm(K*tau). As shown, parameterise a DTMC at lag time tau
T2 = np.eye(n)+pi_arr-np.linalg.inv(np.eye(n)-DT_arr+np.dot(pi_arr,DT_arr)) # alternative expression
for i in range(n): assert abs(np.sum(T1[i,:])-1.)<1.E-08
print "\nT2:\n", T2
