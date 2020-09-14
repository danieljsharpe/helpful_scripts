''' toy python script to demonstrate various computations (e.g. inversion to obtain fundamental matrix for
    both irreducible and absorbing Markov chains), starting from a discrete-time Markov chain and also showing
    the continuous-time case'''

import numpy as np
from math import sqrt
from scipy import linalg

# discrete-time irreducible stochastic transition matrix - the last state is considered to be the target state
P = np.array([[0.5,0.2,0.15,0.15,0.],[0.15,0.75,0.1,0.,0.],[0.2,0.1,0.45,0.1,0.15],
              [0.05,0.,0.25,0.7,0.],[0.,0.,0.1,0.,0.9]],dtype=float)

# committor probabilities when first node is initial state and last node is target state
q = np.array([1.2480468750024063e-01,1.4062499999997560e-01,3.5156250000032622e-01,2.9296875000007163e-01,1.0000000000000000])

'''
P = np.array(
[[  9.54526405e-01,   2.58865114e-04,   5.38032268e-03,   3.91099416e-02,
    0.00000000e+00,   6.77058479e-04,   0.00000000e+00,   0.00000000e+00,
    4.74074113e-05],
 [  1.60395631e-04,   9.62464024e-01,   4.15612113e-03,   0.00000000e+00,
    0.00000000e+00,   4.94312868e-05,   1.55180940e-04,   4.01679762e-03,
    2.89980494e-02],
 [  8.91712862e-04,   1.11169562e-03,   9.80410284e-01,   6.49127552e-05,
    1.73328389e-05,   1.21386689e-02,   5.36169663e-03,   2.51957601e-06,
    1.17704114e-06],
 [  1.35617282e-01,   0.00000000e+00,   1.35812957e-03,   8.35455356e-01,
    0.00000000e+00,   2.75692318e-02,   0.00000000e+00,   0.00000000e+00,
    0.00000000e+00],
 [  0.00000000e+00,   0.00000000e+00,   3.59418846e-04,   0.00000000e+00,
    8.35936792e-01,   5.56892179e-03,   1.57983492e-01,   1.51375900e-04,
    0.00000000e+00],
 [  4.89239686e-04,   5.76472249e-05,   5.29236579e-02,   5.74502989e-03,
    1.17089719e-03,   9.39429229e-01,   1.84298721e-04,   0.00000000e+00,
    0.00000000e+00],
 [  0.00000000e+00,   8.40391489e-05,   1.08554494e-02,   0.00000000e+00,
    1.54250341e-02,   8.55833142e-05,   9.70717218e-01,   2.83267643e-03,
    0.00000000e+00],
 [  0.00000000e+00,   8.32643063e-03,   1.95257931e-05,   0.00000000e+00,
    5.65726875e-05,   0.00000000e+00,   1.08425796e-02,   9.80570039e-01,
    1.84852656e-04],
 [  1.08632034e-04,   1.07241146e-01,   1.62737111e-05,   0.00000000e+00,
    0.00000000e+00,   0.00000000e+00,   0.00000000e+00,   3.29791485e-04,
    8.92304157e-01]])
'''

n = np.shape(P)[0] # number of states (nodes)

for i in range(n): # check row-stochasticity of transition matrix
    assert abs(np.sum(P[i,:])-1.) < 1.E-08
for i in range(n-1): # check first step relation for committor probabilities for non-target nodes
    q_nbr = 0.
    for j in range(n):
        if j==0: continue
        q_nbr += P[i,j]*q[j]
    assert abs(q[i]-q_nbr)<1.E-07

Pr = np.zeros((n,n)) # stochastic matrix for transition process (reactive portion of trajectories)
Pn = np.zeros((n,n)) # for nonreactive process (here, the final node is a dummy node introduced to ensure that the nonreactive traj terminates)
for i in range(n-1):
    for j in range(n):
        if j==0:
            Pr[i,j] = 0.
        else:
            Pr[i,j] = P[i,j]*q[j]/q[i] # here, q_i for i=initial node should satisfy first step relation
        qi=q[i]; qj=q[j]
        if i==0: qi=0.
        if j==0: qj=0.
        if j!=n-1:
            Pn[i,j] = P[i,j]*(1.-qj)/(1.-qi)
        elif i==0:
            Pn[i,j] = q[i]
Pr[-1,-1] = 1.
Pn[-1,-1] = 1.

Q = P[:-1,:-1] # substochastic matrix of nonabsorbing states (i.e. Q-matrix of transition matrix in canonical form)
tau = 0.05 # lag time of DTMC

### INVESTIGATE FUNDAMENTAL MATRIX OF ABSORBING DTMC

''' in the following, note that the element n_{ij} of the fundamental matrix N is equal to the expected number of times
    node j is visited, given that the process starts in node i, prior to absorption. This does not include the fact that
    the process starts in node i, so the row sums are equal to the no. of steps prior to absorption when starting in state i'''

M = np.eye(n-1,dtype=float)-Q # Markovian kernel
N = np.linalg.inv(M) # fundamental matrix of absorbing Markov chain
Nvar = np.dot(N,2.*np.diag(np.diagonal(N))-np.eye(n-1))-(N*N)
H = np.dot(N-np.eye(n-1),np.diag([1./x for x in np.diagonal(N)])) # visitation probabilities along first passage paths
Nr = np.linalg.inv(np.eye(n-1)-Pr[:-1,:-1]) # expected numbers of node visits along transition (reactive) paths
Nn = np.linalg.inv(np.eye(n-1)-Pn[:-1,:-1]) # expected numbers of node visits along nonreactive paths
Eta = np.dot(Nr-np.eye(n-1),np.diag([1./x for x in np.diagonal(Nr)])) # visitation probabilities along transition (reactive) paths
Etn = np.dot(Nn-np.eye(n-1),np.diag([1./x for x in np.diagonal(Nn)])) # visitation probabilities along nonreactive paths

l = np.dot(N,np.ones(n-1)) # vector of expected path lengths
lvar =  np.dot((2.*N)-np.eye(n-1),l)-(l*l)  # vector of variances associated with MFPTs

print "\ndiscrete-time transition probability matrix:\n", P
print "\nstochastic matrix for reactive process:\n", Pr
print "\nstochastic matrix for nonreactive process:\n", Pn
print "\nfundamental matrix of absorbing chain (mean numbers of node visits):\n", N
print "\nvariances in numbers of node visits:\n", Nvar
print "\nvisitation probability matrix:\n", H
print "\nmean numbers of node visits for reactive paths:\n", Nr
print "\nreactive visitation probability matrix:\n", Eta
print "\nmean numbers of node visits for nonreactive paths:\n", Nn
print "\nnonreactive visitation probability matrix:\n", Etn

print "\n"
for i in range(n-1):
    print "expected time to absorption (MFPT) / variance thereof, when starting in state %i:    %.6f   /   %.6f" \
        % (i+1,l[i]*tau,lvar[i]*tau)

### EIGENDECOMPOSITION OF DISCRETE-TIME CHAIN

evals, revecs = np.linalg.eig(P.T) # calculate right eigenvectors of the irreducible stochastic matrix
revecs = np.array([revecs[:,i] for i in evals.argsort()[::-1]]).T
evals, levecs = np.linalg.eig(P) # calculate left eigenvectors of the irreducible stochastic matrix
levecs = np.array([levecs[:,i] for i in evals.argsort()[::-1]]).T
evals = evals[evals.argsort()[::-1]]
print "\neigenvalues\n", evals
#print "\nright eigenvectors:\n", revecs[:,:].T
#print "\nleft eigenvectors:\n", levecs[:,:].T
assert abs(evals[evals.argsort()[-1]]-1.)<1.E-08 # there should be a single dominant eigenvalue equal to unity
pi = revecs[:,0]/np.sum(revecs[:,0]) # equilibrium occupation probabilities (normalised dominant right eigenvector)
# compute normalisation factors
tmp_arr_r = np.zeros((n,n),dtype=float) # diagonal elems are normalisation factors associated with right eigenvectors
tmp_arr_l = np.zeros((n,n),dtype=float) # diagonal elems are normalisation factors associated with left eigenvectors
for i in range(n):
    for j in range(n): 
        for k in range(n):
            tmp_arr_r[i,j] += revecs[k,i]*revecs[k,j]/pi[k]
            tmp_arr_l[i,j] += levecs[k,i]*levecs[k,j]*pi[k]
# normalise
for i in range(n):
    revecs[:,i] *= 1./sqrt(tmp_arr_r[i,i])
    levecs[:,i] *= 1./sqrt(tmp_arr_l[i,i])
#print "\ntmp_arr_r:\n", tmp_arr_r
#print "\ntmp_arr_l:\n", tmp_arr_l
# these are the normalisation checks in Buchete and Hummer J Phys Chem B 2008, for some reason the eigenvectors of the 5-state DTMC do not pass these tests?
# check orthonormality of both left and right eigenvectors
for i in range(n):
#    print "i:", i+1
    assert abs(np.sum(revecs[:,i])-(lambda i: 1. if i==0 else 0.)(i))<1.E-08
    assert abs(abs(np.dot(pi,levecs[:,i]))-(lambda i: 1. if i==0 else 0.)(i))<1.E-08
    for j in range(i,n):
#        print "\t\tdot product of i-th levec with j-th revec:", j+1, "\t", abs(np.dot(revecs[:,i],levecs[:,j]))
#        assert abs(abs(np.dot(revecs[:,i],levecs[:,j]))-(lambda i,j: 1. if i==j else 0.)(i,j))<1.E-08
        if i==j: continue
#        assert abs(tmp_arr_r[i,j])<1.E-07
#        assert abs(tmp_arr_l[i,j])<1.E-07
#        print "\t j:", j+1, "\t", tmp_arr_r[i,j], "       ", tmp_arr_l[i,j]

# compute matrix of all pairwise inter-node MFPTs from eigenspectrum
MFPT_eig = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(n):
        for k in range(1,n):
            MFPT_eig[i,j] += (1./pi[j])*(1.+(evals[k]/(1.-evals[k]))*revecs[j,k]*(levecs[j,k]-levecs[i,k]))
#            MFPT_eig[i,j] += -1.*(tau/(pi[j]*np.log(evals[k])))*revecs[j,k]*(levecs[j,k]-levecs[i,k])

print "\nnormalised right eigenvectors:\n", revecs.T
#print "\nnormalised left eigenvectors:\n", levecs.T
print "\nstationary distribution:\n", pi
print "\nmatrix of MFPTs (from eigenspectrum):\n", MFPT_eig*tau


# INVESTIGATE FUNDAMENTAL MATRIX OF IRREDUCIBLE DTMC

Z = np.linalg.inv(np.eye(n)-P+np.outer(np.ones(n),pi)) # Kemeny and Snell's fundamental matrix. NB may have negative entries
A = Z-np.outer(np.ones(n),pi) # Meyer's group inverse
D = np.diag([1./pi[k] for k in range(n)])

# element-wise computation using fundamental matrix
MFPT_Z = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(n):
        MFPT_Z[i,j] = (Z[j,j]-Z[i,j])
        if i==j: MFPT_Z[i,i] += 1.
        MFPT_Z[i,j] *= 1./pi[j]

MFPT_A = np.dot(np.eye(n)-A+np.dot(np.ones((n,n)),np.diag(np.diagonal(A))),D) # single-line computation using group inverse

for i in range(n-1): assert abs(MFPT_A[i,-1]-l[i])<1.E-07 # check consistency with absorbing formulation for transitions to final node
for mfptz, mfpta in zip(MFPT_Z.flatten(),MFPT_A.flatten()): assert abs(mfptz-mfpta)<1.E-07 # element-wise and single-line methods should give same answer
# check first step relation for MFPTs
MFPT_fsr = np.ones((n,n))
for i in range(n):
    for j in range(n):
        for k in range(n):
            if k==j: continue
            MFPT_fsr[i,j] += P[i,k]*MFPT_A[k,j]
for mfpt_fsr, mfpta in zip(MFPT_fsr.flatten(),MFPT_A.flatten()): assert abs(mfpt_fsr-mfpta)<1.E-07

zeta_P = 1.+np.sum([1./(1.-evals[k]) for k in range(1,n)]) # Kemeny constant from DTMC
zeta_K = np.sum([abs(tau/np.log(evals[k])) for k in range(1,n)]) # Kemeny constant from CTMC
assert abs(zeta_P-np.trace(Z))<1.E-07
#assert abs(zeta_P-zeta_K)<1.E-07
# Kemeny constant can be written as a weighted sum of MFPTs to target nodes (choice of initial node is arbitrary)
for i in range(n): assert abs(np.dot(pi,MFPT_A[i,:])-np.trace(Z))<1.E-07
print "\nthe Kemeny constant is constant for the MFPT matrix computed from eigenspectrum, but value does not match Tr(Z) - eigenvectors not normalised properly?"
for i in range(n): print np.dot(pi,MFPT_eig[i,:])

Gd = np.diag(np.diagonal(np.dot(np.eye(n)-np.outer(np.ones(n),pi),A)))
Mvd = D+(2.*np.dot(np.dot(D,Gd),D))
Var_A = 2.*(np.dot(A,MFPT_A)-np.dot(np.ones(n),np.diag(np.diagonal(np.dot(A,MFPT_A)))))
Var_A += np.dot(np.dot(MFPT_A,np.diag(pi)),Mvd)
Var_A = Var_A-(MFPT_A*MFPT_A) # variances in FPT distributions for transitions between all pairs of nodes (i,j)
for i in range(n-1): assert abs(Var_A[i,-1]-lvar[i])<1.E-07 # check consistency with absorbing formulation for transitions to final node

print "\ngroup inverse (aka deviation) matrix:\n", A
print "\nKemeny constant:", np.trace(A)+1.
print "\nmatrix of MFPTs (from group inverse):\n", MFPT_A*tau
print "\nvariances of FPT distribution (from group inverse):\n", Var_A*tau

# COMPUTE CONTINUOUS-TIME MARKOV CHAIN AND ANALOGOUS QUANTITIES TO THOSE GIVEN ABOVE

MFPT_A_ctmc = (MFPT_A-np.diag(np.diagonal(MFPT_A)))*tau # MFPT matrix in continuous-time
print "\nMFPT_A_ctmc:\n", MFPT_A_ctmc
K = np.dot(np.linalg.inv(MFPT_A_ctmc),D-np.ones((n,n))) # transition rate matrix

#pi_arr = np.outer(pi,np.ones(n))
#K = pi_arr-np.linalg.inv(pi_arr-np.dot(np.dot(np.diag(pi),MFPT_A_ctmc),np.eye(n)-pi_arr))

print "\n"
for i in range(n): print "rate matrix sum of row:", i, np.sum(K[i,:])

print "\ntransition rate matrix:\n", K
T_rew = linalg.expm(K*tau)
print "\nrecovered transition probability matrix:\n", T_rew

#pi_mfpt = np.dot(np.diag(pi),MFPT_A)
#T_new = np.dot(np.linalg.inv(np.eye(n)-pi_mfpt),pi_arr-pi_mfpt)
#print "\nrecovered stochastic matrix:\n", T_new

'''
print "\ntest:\n", (np.dot(MFPT_A,pi)-1.)*tau
print "\nshould be identical to:\n", (np.dot(P,np.dot(MFPT_A,pi))-1.)*tau
print "\ntest 2:\n", np.dot(MFPT_A_ctmc,pi)
'''
