''' toy python script to demonstrate various computations (e.g. inversion to obtain fundamental matrix for
    both irreducible and absorbing Markov chains), starting from a discrete-time Markov chain and also showing
    the continuous-time case'''

from __future__ import print_function
import numpy as np
from math import sqrt
from scipy import linalg

''' basic implementation of Daniel's state reduction algorithm to compute the fundamental matrix for a reducible
    Markov chain. Args: (stochastic matrix, no. of transient states) '''
def stateredn_absfundmtx(T,nQ):
    nS = np.shape(T)[0] # total no. of states
#    print("T shape:",np.shape(T))
    # construct augmented matrix
    Naug = np.pad(T,(0,nQ),mode="constant",constant_values=((0.)))
    ''' NB zero weights for transitions from absorbing to dummy nodes and vice versa, and for transitions between
        dummy nodes, but weight=unity for transitions to and from a dummy node and its transient partner '''
    for i in range(nQ):
        Naug[i,i+nS], Naug[i+nS,i] = 1., 1.
#    print("Naug shape:",np.shape(Naug))
#    print("\ninitial augmented matrix:\n",Naug)
    for x in range(nQ): # eliminate each of the transient nodes
        Naug[x,x]=0.
        # calculate GT factor in numerically stable manner
        Nfac = 0.
        for j in range(x+1,nS):
            Nfac += Naug[x,j]
        for i in range(x+1,nS+nQ):
            for j in range(x+1,nS+nQ):
                Naug[i,j] += Naug[i,x]*Naug[x,j]/Nfac
        Naug[x,:], Naug[:,x] = 0., 0. # remove row and column corresponding to eliminated node
    return Naug[nS:,nS:]


### SETUP AND BASIC CHECKS OF INITIAL DISCRETE-TIME MARKOV CHAIN ###

n = 15 # no. of states (nodes)
nB = 6 # no. of nodes in initial state (must be listed as first nodes)
ndB = 3 # no. of nodes at the boundary of the initial state B (must be listed after the internal nodes of the initial state
nA = 3 # no. of nodes in absorbing state (must be listed as the final nodes)
do_ctmc = False # if True, DTMC from first code section is transformed to a CTMC and further operations are performed
use_statereduction = True # fundamental matrix of absorbing MC is calcd by state reduction (T) or by simple matrix inversion operation (F)

'''
# discrete-time irreducible stochastic transition matrix - the last state is considered to be the target state
# note: this stochastic matrix satisfies global, but not detailed, balance
P = np.array([[0.50,0.20,0.15,0.15,0.00],
              [0.15,0.75,0.10,0.00,0.00],
              [0.20,0.10,0.45,0.10,0.15],
              [0.05,0.00,0.25,0.70,0.00],
              [0.00,0.00,0.10,0.00,0.90]],dtype=float)
tau = 0.05 # lag time of DTMC
tau_vec = np.array([tau]*np.shape(P)[0])
# committor probabilities when first node is initial state and last node is target state
# NB for qp vector, committor probs for *all* (incl initial) nodes must satisfy first-step relation
#    for q  vector, committor probs for initial nodes are zero by definition (initialised below)
qp = np.array([1.2480468750024063e-01,1.4062499999997560e-01,3.5156250000032622e-01,2.9296875000007163e-01,1.0000000000000000])
'''

#'''
# alternatively, read in data from pickled branching probability matrix, vector of mean waiting times,
#and A<-B committor probabilities. (NB: the computed variances in the FPT distribution are then not valid)
P = np.load("branchmtx.pkl")
tau_vec = np.load("meanwaittimes.pkl")
qp = np.load("committors.pkl")
#'''

#'''
# quack - need to do some reordering of initial nodes for 15-state CTMC model
P[[3,4]] = P[[4,3]] # swap row
P[:,[3,4]] = P[:,[4,3]] # swap columns
P[[1,3]] = P[[3,1]]
P[:,[1,3]] = P[:,[3,1]]
tau_vec[3], tau_vec[4] = tau_vec[4], tau_vec[3]
tau_vec[1], tau_vec[3] = tau_vec[3], tau_vec[1]
qp[3], qp[4] = qp[4], qp[3]
qp[1], qp[3] = qp[3], qp[1]
#'''

assert n==np.shape(P)[0]
for i in range(n): # check row-stochasticity of transition matrix
    assert abs(np.sum(P[i,:])-1.)<1.E-08

# check first step relation for committor probabilities for non-target nodes
for i in range(n-nA):
    q_nbr = 0.
    for j in range(n):
        if j<nB: continue # ignore transitions to nodes of initial set
        q_nbr += P[i,j]*qp[j]
    assert abs(qp[i]-q_nbr)<1.E-08
q = qp.copy()
for i in range(nB): q[i]=0.


### INVESTIGATE FUNDAMENTAL MATRIX OF THE ABSORBING (REDUCIBLE) MARKOV CHAIN

Pr = np.zeros((n,n)) # stochastic matrix for transition process (reactive portion of trajectories)
Pn = np.zeros((n,n)) # for nonreactive process (here, the final node is a dummy node introduced to ensure that the nonreactive traj terminates)
for i in range(n-nA):
    for j in range(n):
        # reactive (transition) process
        if j<nB: # reactive process does not re-enter initial state
            Pr[i,j] = 0.
        elif i>=(nB-ndB): # only initial nodes that are at the boundary of the B state are relevant to the reactive process
            Pr[i,j] = P[i,j]*qp[j]/qp[i] # here, committor prob for i=initial node should satisfy first step relation
        # nonreactive process
        if j<n-nA: # j is not an absorbing state
            Pn[i,j] = P[i,j]*(1.-q[j])/(1.-q[i]) # here, committor prob should be 0 for initial states, *not* satisfying first step relation
        elif i<nB: # i is an initial state
            Pn[i,j] = qp[i]
for i in range(n-nA,n): # only outgoing transitions for absorbing states are self-loops
    Pr[i,i] = 1.
    Pn[i,i] = 1.

Q = P[:-nA,:-nA] # substochastic matrix of nonabsorbing states (i.e. Q-matrix of transition matrix in canonical form)

''' in the following, note that the element n_{ij} of the fundamental matrix N is equal to the expected number of times
    node j is visited, given that the process starts in node i, prior to absorption. This does not include the fact that
    the process starts in node i, so the row sums are equal to the no. of steps prior to absorption when starting in state i'''

M = np.eye(n-nA,dtype=float)-Q # Markovian kernel
if not use_statereduction: # fundamental matrices of absorbing Markov chain computed by simple matrix inversion operations
    N = np.linalg.inv(M)
    Nr = np.linalg.inv(np.eye(n-nB+ndB-nA)-Pr[nB-ndB:-nA,nB-ndB:-nA]) # expected numbers of node visits along transition (reactive) paths
    Nn = np.linalg.inv(np.eye(n-nA)-Pn[:-nA,:-nA]) # expected numbers of node visits along nonreactive paths
else: # fundamental matrices of absorbing Markov chain computed by state reduction algorithm
    N = stateredn_absfundmtx(P,n-nA)
    Nr = stateredn_absfundmtx(Pr[nB-ndB:,nB-ndB:],n-nB+ndB-nA)
#    print("Nr shape:",np.shape(Nr))
#    quit()
    Nn = stateredn_absfundmtx(Pn,n-nA)
Nvar = np.dot(N,2.*np.diag(np.diagonal(N))-np.eye(n-nA))-(N*N)
H = np.dot(N-np.eye(n-nA),np.diag([1./x for x in np.diagonal(N)])) # visitation probabilities along first passage paths
Hr = np.dot(Nr-np.eye(n-nB+ndB-nA),np.diag([1./x for x in np.diagonal(Nr)])) # visitation probabilities along transition (reactive) paths
Hn = np.dot(Nn-np.eye(n-nA),np.diag([1./x for x in np.diagonal(Nn)])) # visitation probabilities along nonreactive paths
# pad matrices corresponding to reactive quantities, if there are internal (i.e. non-boundary) initial nodes
if nB!=ndB:
    Nr = np.pad(Nr,(nB-ndB,0),mode="constant")
    Hr = np.pad(Hr,(nB-ndB,0),mode="constant")

l = np.dot(N,np.ones(n-nA)) # vector of expected first passage path lengths
m = np.dot(N,tau_vec[:n-nA]) # vector of mean first passage times (MFPTs)
lvar =  np.dot((2.*N)-np.eye(n-nA),l)-(l*l)  # vector of variances in first passage path lengths
mvar = lvar*tau_vec[:n-nA]*tau_vec[:n-nA] # vector of variances associated with MFPTs (only valid for DTMC)

print("\ntransition probability matrix:\n",P)
print("\nlag time (if DTMC) or mean waiting time (if CTMC) vector:\n",tau_vec)
print("\ncommittor probability vector:\n",qp)
print("\nstochastic matrix for reactive process:\n",Pr)
print("\nstochastic matrix for nonreactive process:\n",Pn)
print("\nfundamental matrix of absorbing chain (mean numbers of node visits):\n",N)
print("\nvariances in numbers of node visits:\n",Nvar)
print("\nvisitation probability matrix:\n",H)
print("\nmean numbers of node visits for reactive paths:\n",Nr)
print("\nreactive visitation probability matrix:\n",Hr)
print("\nmean numbers of node visits for nonreactive paths:\n",Nn)
print("\nnonreactive visitation probability matrix:\n",Hn)

print("\n\nnode  /     MFPT      /   FPT variance")
for i in range(n-nA): print("{:4d}".format(i+1),"\t","{:.6e}".format(m[i]),"\t","{:.6e}".format(mvar[i]))


### EIGENDECOMPOSITION OF DISCRETE-TIME MARKOV CHAIN

evals, revecs = np.linalg.eig(P.T) # calculate right eigenvectors of the irreducible stochastic matrix
revecs = np.array([revecs[:,i] for i in evals.argsort()[::-1]])
evals, levecs = np.linalg.eig(P) # calculate left eigenvectors of the irreducible stochastic matrix
levecs = np.array([levecs[:,i] for i in evals.argsort()[::-1]])
evals = evals[evals.argsort()[::-1]]
assert abs(evals[evals.argsort()[-1]]-1.)<1.E-08 # there should be a single dominant eigenvalue equal to unity
pi = np.real(revecs[0,:]/np.sum(revecs[0,:])) # equilibrium occupation probabilities (normalised dominant right eigenvector)
print("\nstationary distribution:\n",pi)
print("\neigenvalues:\n",evals)


### DANIEL TRANSITION PATH THEORY QUANTITIES

# initial probability distribution for first passage trajectories is local equilibrium distribution within B
p0tot = 0.
for i in range(nB):
    p0tot += pi[i]
p0 = pi[:nB].copy()/p0tot
# expected number of times that transient nodes are visited along A<-B first passage trajectories
theta = np.dot(p0,N[:nB,:])
# expected number of times that transient nodes are visited along nonreactive B<-B trajectories
thetanr = np.zeros(n-nA,dtype=float)
for i in range(n-nA):
    thetanr[i] = theta[i]*(1.-q[i]) # committor probability must be zero for initial state in this context
# initial probability distribution for *reactive* trajectories (contained at boundary of B state)
p0r = np.zeros(nB,dtype=float)
for i in range(nB):
    for j in range(n):
        p0r[i] += P[i,j]*q[j]  # committor probability must be zero for initial state in this context
    p0r[i] *= thetanr[i]
assert abs(np.sum(p0r)-1.)<1.E-08
# expected number of times that transient nodes are visited along reactive A<-B (i.e. transition) trajectories
thetar = np.dot(p0r,Nr[:nB,:])
# initial probability distribution for *reactive* trajectories at steady-state
p0rss = np.zeros(nB,dtype=float)
for i in range(nB):
    for j in range(n):
        p0rss[i] += P[i,j]*q[j] # committor probability must be zero for initial state in this context
    p0rss[i] *= pi[i]
p0rss *= 1./np.sum(p0rss)
# expected number of times that transient nodes are visited along reactive trajectories at steady-state
thetarss = np.dot(p0rss,Nr[:nB,:])
# visitation probabilities for transient nodes
eta = np.dot(p0,H[:nB,:]) # visitation probabilities for first passage trajectories
etar = np.dot(p0r,Hr[:nB,:]) # visitation probabilities for reactive trajectories
etarss = np.dot(p0rss,Hr[:nB,:]) # visitation probabilities for reactive trajectories at steady-state

print("\ninitial probability distribution for reactive trajectories:\n",p0r)
print("initial probability distribution for reactive trajectories at steady-state:\n",p0rss)
print("\n\nexpected number of node visits for transient states:\n", \
      "node /  first passage  /   reactive   /  nonreactive  /  reactive steady-state")
for i in range(n-nA):
    print("{:4d}".format(i+1),"\t","{:.6e}".format(theta[i]),"\t","{:.6e}".format(thetar[i]),"\t", \
          "{:.6e}".format(thetanr[i]),"\t","{:.6e}".format(thetarss[i]))
print("\n\nvisitation probabilities for transient states:\n", \
      "node /  first passage  /   reactive   /  reactive steady-state")
for i in range(n-nA):
    print("{:4d}".format(i+1),"\t","{:.6e}".format(eta[i]),"\t","{:.6e}".format(etar[i]),"\t", \
          "{:.6e}".format(etarss[i]))
print("\nA<-B MFPT:\t","{:.6e}".format(np.dot(p0,m[:nB])))

quit()


# INVESTIGATE FUNDAMENTAL MATRIX OF IRREDUCIBLE DTMC

Zinv = np.eye(n)-P+np.outer(np.ones(n),pi)
Z = np.linalg.inv(Zinv) # Kemeny and Snell's fundamental matrix. NB may have negative entries
A = Z-np.outer(np.ones(n),pi) # Meyer's group inverse
D = np.diag([1./pi[k] for k in range(n)])
print("\ncondition number in matrix inversion to obtain Z:",np.linalg.cond(Zinv))

# element-wise computation using fundamental matrix
MFPT_Z = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(n):
        MFPT_Z[i,j] = Z[j,j]-Z[i,j]
        if i==j: MFPT_Z[i,i] += 1.
        MFPT_Z[i,j] *= 1./pi[j]

MFPT_A = np.dot(np.eye(n)-A+np.dot(np.ones((n,n)),np.diag(np.diagonal(A))),D) # single-line computation using group inverse

for i in range(n-1): assert abs(MFPT_A[i,-1]-l[i])<1.E-07 # check consistency with absorbing formulation for transitions to final node
for mfptz, mfpta in zip(MFPT_Z.flatten(),MFPT_A.flatten()): assert abs(mfptz-mfpta)<1.E-08 # element-wise and single-line methods should give same answer
# check first step relation for MFPTs
MFPT_fsr = np.ones((n,n))
for i in range(n):
    for j in range(n):
        for k in range(n):
            if k==j: continue
            MFPT_fsr[i,j] += P[i,k]*MFPT_A[k,j]
for mfpt_fsr, mfpta in zip(MFPT_fsr.flatten(),MFPT_A.flatten()): assert abs(mfpt_fsr-mfpta)<1.E-08

zeta_P = 1.+np.sum([1./(1.-evals[k]) for k in range(1,n)]) # Kemeny constant from DTMC
zeta_K = np.sum([abs(tau/np.log(evals[k])) for k in range(1,n)]) # Kemeny constant from CTMC
assert abs(zeta_P-np.trace(Z))<1.E-08
#assert abs(zeta_P-zeta_K)<1.E-08
# Kemeny constant can be written as a weighted sum of MFPTs to target nodes (choice of initial node is arbitrary)
for i in range(n): assert abs(np.dot(pi,MFPT_A[i,:])-np.trace(Z))<1.E-08

Gd = np.diag(np.diagonal(np.dot(np.eye(n)-np.outer(np.ones(n),pi),A)))
Mvd = D+(2.*np.dot(np.dot(D,Gd),D))
Var_A = 2.*(np.dot(A,MFPT_A)-np.dot(np.ones(n),np.diag(np.diagonal(np.dot(A,MFPT_A)))))
Var_A += np.dot(np.dot(MFPT_A,np.diag(pi)),Mvd)
Var_A = Var_A-(MFPT_A*MFPT_A) # variances in FPT distributions for transitions between all pairs of nodes (i,j)
for i in range(n-1): assert abs(Var_A[i,-1]-lvar[i])<2.E-05 # check consistency with absorbing formulation for transitions to final node

print("\ngroup inverse (aka deviation) matrix:\n",A)
print("\nKemeny constant:",np.trace(A)*tau)
print("\nmatrix of MFPTs (from group inverse):\n",MFPT_A*tau)
print("\nvariances of FPT distribution (from group inverse):\n",Var_A*(tau**2))

if not do_ctmc: quit()


# COMPUTE CONTINUOUS-TIME MARKOV CHAIN AND ANALOGOUS QUANTITIES TO THOSE GIVEN ABOVE

MFPT_A_ctmc = (MFPT_A-np.diag(np.diagonal(MFPT_A)))*tau # MFPT matrix in continuous-time
K = np.dot(D-np.ones((n,n)),np.linalg.inv(MFPT_A_ctmc)) # transition rate matrix
for i in range(n): assert abs(np.sum(K[i,:]))<1.E-08
evals_K, revecs_K = np.linalg.eig(K.T)
revecs_K = np.array([revecs_K[:,i] for i in evals_K.argsort()[::-1]])
pi_K = revecs_K[0,:]/np.sum(revecs_K[0,:])
for i in range(n): assert abs(pi[i]-pi_K[i])<1.E-08
T_rew = linalg.expm(K*tau)
tau_vec = np.array([1./abs(K[i,i]) for i in range(n)]) # vector of mean waiting times
B = np.zeros((n,n),dtype=float) # branching probability matrix
for i in range(n):
    for j in range(i+1,n):
        B[i,j] = K[i,j]/abs(K[i,i])
        B[j,i] = K[j,i]/abs(K[j,j])

print("\ntransition rate matrix:\n",K)
print("\nrecovered transition probability matrix:\n",T_rew) # NB is not necessary same as transition probability matrix for original DTMC

NK = np.linalg.inv(np.eye((n-1),dtype=float)-B[:-1,:-1]) # fundamental matrix of absorbing CTMC
MFPT_NK = np.dot(NK,tau_vec[:-1]) # MFPTs for transitions to final state in CTMC
NKvar = np.dot(NK,2.*np.diag(np.diagonal(NK))-np.eye(n-1))-(NK*NK) # variances in number of node visits
#tau_var = np.tile(tau_vec[:-1]*tau_vec[:-1],n-1).reshape((n-1,n-1)) # variance of exp distribn is square of mean
#Xvar = (NK*NK*tau_var)+(tau_var*NKvar)+(tau_var*NKvar)
#Var_NK = np.dot(Xvar,np.ones(n-1))
Var_NK = np.dot(NKvar,(tau_vec[:-1]*tau_vec[:-1]))

ZK = np.linalg.inv(np.outer(np.ones(n),pi)-K)-np.outer(np.ones(n),pi) # fundamental matrix of irreducible CTMC
MFPT_ZK = np.dot(np.eye(n)-ZK+np.dot(np.ones((n,n)),np.diag(np.diagonal(ZK))),D)
for i in range(n-1): assert abs(MFPT_NK[i]-MFPT_ZK[i,-1])<1.E-08

print("\nMFPTs for CTMC:\n",MFPT_ZK) # MFPTs for CTMC should be the same as for DTMC from which CTMC was derived
print("\nNK:\n",NK,"\nNKvar:\n",NKvar)
#print("\nXvar:\n",Xvar)
print("\nvariances of FPT distribution for CTMC:\n",Var_NK)
print("\nKemeny constant for CTMC:\n",np.trace(ZK))

# check if Markov chain satisfies the *detailed* balance condition
reversible=True
for i in range(n):
    for j in range(i+1,n):
        if abs((P[i,j]*pi[i])-(P[j,i]*pi[j]))>1.E-08: reversible=False
if not reversible: quit()

### FOR A REVERSIBLE MARKOV CHAIN, COMPUTE NORMALIZED EIGENVECTORS AND COMPUTE MFPT FROM EIGENSPECTRUM

tmp_arr_r = np.zeros((n,n),dtype=float) # diagonal elems are normalisation factors associated with right eigenvectors
tmp_arr_l = np.zeros((n,n),dtype=float) # diagonal elems are normalisation factors associated with left eigenvectors
for i in range(n):
    for j in range(n):
        for k in range(n):
            tmp_arr_r[i,j] += revecs[i,k]*revecs[j,k]/pi[k]
            tmp_arr_l[i,j] += levecs[i,k]*levecs[j,k]*pi[k]
# normalise
for i in range(n):
    revecs[i,:] *= 1./sqrt(tmp_arr_r[i,i])
    levecs[i,:] *= 1./sqrt(tmp_arr_l[i,i])

# these are the orthonormality conditions in Buchete and Hummer J Phys Chem B 2008, some of which only apply if *detailed* balance is satisfied
for i in range(n):
#    print("i:",i+1)
    assert abs(np.sum(revecs[i,:])-(lambda i: 1. if i==0 else 0.)(i))<1.E-08
    assert abs(abs(np.dot(pi,levecs[i,:]))-(lambda i: 1. if i==0 else 0.)(i))<1.E-08
    for j in range(i,n):
#        print("\tdot product of i-th levec with j-th revec:", j+1, "\t", abs(np.dot(revecs[i,:],levecs[j,:])))
        assert abs(abs(np.dot(revecs[i,:],levecs[j,:]))-(lambda i,j: 1. if i==j else 0.)(i,j))<1.E-08
        if not reversible or i==j: continue
        assert abs(tmp_arr_r[i,j])<1.E-08
        assert abs(tmp_arr_l[i,j])<1.E-08
#        print("\t j:", j+1, "\t", tmp_arr_r[i,j], "       ", tmp_arr_l[i,j])

# compute matrix of all pairwise inter-node MFPTs from eigenspectrum
MFPT_eig = np.zeros((n,n),dtype=float)
for i in range(n):
    for j in range(n):
        for k in range(1,n):
            MFPT_eig[i,j] += (1./pi[j])*(evals[k]/(1.-evals[k]))*revecs[j,k]*(levecs[j,k]-levecs[i,k])
#            MFPT_eig[i,j] += -1.*(tau/(pi[j]*np.log(evals[k])))*revecs[j,k]*(levecs[j,k]-levecs[i,k])
        MFPT_eig[i,j] += 1./pi[j]
for i in range(n): assert abs(np.dot(pi,MFPT_eig[i,:])-np.trace(Z))<1.E-08

print("\nnormalized right eigenvectors:\n",revecs)
print("\nnormalized left eigenvectors:\n",levecs)
print("\nmatrix of MFPTs (from eigenspectrum):\n",MFPT_eig*tau)
