'''
Playing with linear algebra in Python
'''

import numpy as np
from scipy.sparse.linalg import eigs as eigs_iram
import scipy.linalg

### SET PARAMS
n_nodes = 998
n_edges = 3981 # no. of bidirectional edges
n_eigs = 8
d=0.5


ts_wts = np.zeros(2*n_edges,dtype=np.float128)
k_mtx = np.zeros((n_nodes,n_nodes),dtype=np.float128)

with open("ts_weights.dat","r") as ts_wts_f:
    lines = [line.split()[0] for line in ts_wts_f.readlines()]
    for i in range(0,len(lines),2):
        ts_wts[i] = np.exp(np.float128(lines[i]))
        ts_wts[i+1] = np.exp(np.float128(lines[i+1]))

with open("ts_conns.dat","r") as ts_conns_f:
    for i, line in enumerate(ts_conns_f.readlines()):
        line = line.split()
        if int(line[0])==int(line[1]): continue # this indicates a "dead" TS
#        print 2*i, int(line[0]), int(line[1]), np.log(ts_wts[2*i])
#        print (2*i)+1, int(line[1]), int(line[0]), np.log(ts_wts[(2*i)+1])
        k_mtx[int(line[0])-1,int(line[1])-1] = ts_wts[2*i]
        k_mtx[int(line[1])-1,int(line[0])-1] = ts_wts[(2*i)+1]

for i in range(n_nodes):
    k_mtx[i,i] = -np.sum(k_mtx[i,:])
    assert abs(np.sum(k_mtx[i,:]))<1.E-10
    k_mtx[i,i] -= d

'''
print "eigenvalue calculation, standard method, dense"
evecs, evals = np.linalg.eig(k_mtx)
evals = np.sort(evals)
'''

'''
# implicitly restarted Arnoldi method
evals, evecs = eigs_iram(k_mtx.astype(np.float64),n_eigs,which="SM")
evals = np.sort(evals)
print "evals:\n", evals
'''

#'''
# scipy linalg dense method. Can use np.float128
evals, evecs = scipy.linalg.eig(k_mtx)
evals = np.sort(evals)
print "evals:\n", evals
print "length:", np.shape(evals)
#'''
