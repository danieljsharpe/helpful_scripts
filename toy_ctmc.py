import numpy as np
from scipy import linalg
from math import sqrt

'''
# note: this rate matrix satisfies global, but not detailed, balance
K = np.array([[-5.,2.0,2.0,1.0],
              [1.5,-6.,1.5,3.0],
              [3.0,4.0,-9.,2.0],
              [0.5,0.5,1.0,-2.]])
tau = 0.05 # lag time
'''

#'''
# note: this rate matrix satisfies detailed balance
K = np.array([[-0.142479765469,  0.082084998624,    0.030197383422,  0.030197383422,  0.            ],
              [ 0.367879441171,  -0.59100960132,    0.223130160148,  0.,              0.            ],
              [ 0.082084998624,  0.135335283237,   -0.225593824737,  0.004086771438,  0.004086771438],
              [ 0.367879441171,  0.            ,    0.018315638889, -0.38619508006,   0.            ],
              [ 0.,              0.,                0.082084998624,  0.,             -0.082084998624]])
tau = 5. # lag time
#'''

n = np.shape(K)[0]

evals_K, revecs_K = np.linalg.eig(K.T)
idx = evals_K.argsort()[-1]
pi = revecs_K[:,idx]/np.sum(revecs_K[:,idx])
evals_K = evals_K[evals_K.argsort()[::-1]]
D = np.diag([1./pi[i] for i in range(n)])

assert abs(np.sum(pi)-1.)<1.E-08
for i in range(n): assert abs(np.dot(pi,K)[i])<1.E-08 # check global balance condition
print "\ntransition rate matrix:\n", K
print "\nstationary distribution:\n", pi
print "\neigenvalues:\n", evals_K

reversible=True
for i in range(n):
    for j in range(n):
        if abs((K[i,j]*pi[i])-(K[j,i]*pi[j]))>1.E-08: reversible=False
print "\ndetailed balance condition satisifed?", reversible

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

T = linalg.expm(K*tau) # discrete-time transition probability matrix
evals_T, revecs_T = np.linalg.eig(T.T)
idx = evals_T.argsort()[-1]
piT = revecs_T[:,idx]/np.sum(revecs_T[:,idx])
evals_T = evals_T[evals_T.argsort()[::-1]]
for i in range(n): assert abs(piT[i]-pi[i])<1.E-08 # check stationary distribution is identical to CTMC result
print "\ntransition probability matrix:\n", T

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
'''
pi_arr = np.outer(np.ones(n),pi)
DT_arr = np.dot((MFPT/tau)+D,np.diag(pi)) # could instead use MFPT_T to get back expm(K*tau). As shown, parameterise a DTMC at lag time tau
T2 = np.eye(n)+pi_arr-np.linalg.inv(np.eye(n)-DT_arr+np.dot(pi_arr,DT_arr)) # alternative expression
for i in range(n): assert abs(np.sum(T2[i,:])-1.)<1.E-08
print "\nT2:\n", T2
'''

dS = np.zeros((n,n),dtype=float) # matrix of entropy flow along edges (as multiple of Boltzmann constant)
for i in range(n):
    for j in range(i+1,n):
        if K[i,j]==0. and K[j,i]==0.: continue # nodes are not connected in either direction
        dS[i,j] = -np.log(K[i,j]/K[j,i])
        dS[j,i] = -np.log(K[j,i]/K[i,j])
print "\nentropy flow matrix:\n", dS

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
print "\neigenvectors of symmetrised rate matrix:\n", evecs_sym
