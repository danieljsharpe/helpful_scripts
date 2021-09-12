'''
Assuming that all nodes listed in nodes.A belong to one community, and all other nodes to a second community, write a communities.dat file

Daniel J. Sharpe
Feb 2020
'''

import sys

n_nodes = int(sys.argv[1])

nodeA_nodes_flag = [0]*n_nodes
with open("nodes.A","r") as nodea_f:
    for line in nodea_f.readlines():
        nodeA_nodes_flag[int(line.split()[0])-1]=1

with open("communities.dat","w") as comms_f:
    for i in range(n_nodes):
        comms_f.write("%i\n" % nodeA_nodes_flag[i])
