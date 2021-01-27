'''
Read "communities.dat" file and merge two desired communities, and update new community IDs to be consistent
'''

import sys

comm1 = int(sys.argv[1])
comm2 = int(sys.argv[2])

if comm1>comm2:
    comm_tmp = comm2
    comm2 = comm1
    comm1 = comm_tmp

new_comms = []

with open("communities.dat","r") as comms_f:
    for line in comms_f.readlines():
        comm_id = int(line.split()[0])
        if comm_id==comm2: comm_id = comm1
        elif comm_id>comm2: comm_id-=1
        new_comms.append(comm_id)
with open("communities_merged.dat","w") as commsmg_f:
    for comm_id in new_comms:
        commsmg_f.write("%i\n" % comm_id)
