'''
Assuming that all nodes listed in min.A belong to one community, and all other nodes to a second community, write a communities.dat file

Daniel J. Sharpe
Feb 2020
'''

import sys

n_nodes = int(sys.argv[1])

minA_nodes_flag = [0]*n_nodes
with open("min.A","r") as mina_f:
    for line in mina_f.readlines():
        minA_nodes_flag[int(line.split()[0])-1]=1

with open("communities.dat","w") as comms_f:
    for i in range(n_nodes):
        comms_f.write("%i\n" % minA_nodes_flag[i])
