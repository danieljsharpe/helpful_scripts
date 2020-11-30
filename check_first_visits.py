'''
Python script to detect when the first node of a network is visited more than once
Reads from visits.k.dat files, containing numbers of node visits for k-th shortest path
'''

from __future__ import print_function

npaths=1000 # number of paths
for k in range(npaths):
    with open("visits."+str(k+1)+".dat") as visits_f:
        x = int(visits_f.readline()) # number of times first node is visited along path
        if x>1:
            print("first node visited %i times in %i-th shortest path" %(x,k))
