'''
Python script to parse a "minima_groups" file produced by REGROUPFREE in PATHSAMPLE, and write a Daniel-style communities.dat file
and an attractors.dat file containing represenative nodes of each community

Daniel J. Sharpe
'''

from collections import OrderedDict
import sys

mingps_fname = sys.argv[1]
n_nodes = int(sys.argv[2])

attractors = OrderedDict()
comm_ids = [-1]*n_nodes
curr_comm_id = 0
with open(mingps_fname,"r") as mingps_f:
    for line in mingps_f.readlines():
        line = line.split()
        if not line: continue
        if line[0]=="group":
            curr_comm_id += 1
            continue
        comm_ids[int(line[0])-1] = curr_comm_id
        attractors[curr_comm_id] = int(line[0])

with open("communities_regroup.dat","w") as comms_f:
    for comm_id in comm_ids:
        assert comm_id>=0
        comms_f.write("%i\n" % comm_id)
with open("attractors_regroup.dat","w") as atts_f:
    for key, item in attractors.items():
        atts_f.write("%i\n" % item)
