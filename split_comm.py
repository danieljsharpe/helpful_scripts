'''
Read file "communities.dat" and a file "min.comm" of nodes to retain their community ID. All other nodes in the specified community are put in
a new singleton community
'''

import sys

ncomms=int(sys.argv[1])
nnodes=int(sys.argv[2])
activecomm=int(sys.argv[3])

comm_ids=[-1]*nnodes
with open("communities.dat","r") as comms_f:
    for i, line in enumerate(comms_f.readlines()):
        comm_ids[i]=int(line.split()[0])
nodemask=[False]*nnodes
with open("min.comm","r") as min_f:
    for line in min_f.readlines():
        node_id=int(line.split()[0])
        assert comm_ids[node_id-1]==activecomm
        nodemask[node_id-1]=True
for i in range(nnodes):
    if comm_ids[i]!=activecomm or nodemask[i]: # comm ID of node does not change
        continue
    else:
        comm_ids[i]=ncomms
        ncomms+=1
with open("communities_split.dat","w") as commsnew_f:
    for comm_id in comm_ids:
        commsnew_f.write("%i\n" % comm_id)
print "the number of communities is now: %i" % ncomms
