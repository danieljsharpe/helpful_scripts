'''
Read a communities.dat file and get representative (dummy "attractor") nodes for each community
'''

import sys

n_comms = int(sys.argv[1])

attractors = [-1]*n_comms
with open("communities.dat","r") as comms_f:
    for i, line in enumerate(comms_f.readlines()):
        if attractors[int(line.split()[0])]==-1:
            attractors[int(line.split()[0])] = i+1
with open("attractors_dummy.dat","w") as atts_f:
    for att in attractors:
        assert att>0
        atts_f.write("%i\n" % att)
