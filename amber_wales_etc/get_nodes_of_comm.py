'''
Python script to write a file "min.comm" of nodes belonging to specified community
Communities are read from Daniel-style file "communities.dat"

Daniel J. Sharpe
'''

import sys

comm_id=int(sys.argv[1])

nodelist=[]
with open("communities.dat","r") as comms_f:
    for i, line in enumerate(comms_f.readlines()): # NB comms_f is single-col and no header
        if int(line)==comm_id: nodelist.append(i+1)

with open("min.comm","w") as nodes_f:
    for node_id in nodelist:
        nodes_f.write("%i\n" % node_id)
print "written %i nodes to file" % len(nodelist)
