'''
Script to parse a "ts.data" file and (in-place) delete duplicate transitions state (i.e. TSs connecting
the same pairs of minima). The lowest-energy TS connecting any given pair of minima is always retained.

Usage:
python remove_dup_ts.py <ts.data file>
'''

import subprocess
import sys

tsd_fname = sys.argv[1]

# dict-of-dicts representation of KTN, format: G[v1] : {v2: ts_id, ...}, where v1 < v2
# G is used to check if a pair of minima are already connected
G = {}
ts_conns = []
ts_energies = []

# read in data
nts = 0
with open(tsd_fname,"r") as tsd_f:
    for line in tsd_f.readlines():
        line = line.split()
        ts_energies.append(float(line[0]))
        ts_conns.append([int(line[3]),int(line[4])])
        nts += 1

print "Number of transition states: %i" % nts
raw_ts_ids = [i for i in range(1,nts+1)] # used to keep track of line numbers
# construct data structure for KTN, checking for duplicate TSs
for i, ts_conn in enumerate(ts_conns):
    if ts_conn[0] == ts_conn[1]: continue # dead TS
    if ts_conn[0] < ts_conn[1]: # find the lower- and higher-index nodes
        v1, v2 = ts_conn[0], ts_conn[1]
    else:
        v1, v2 = ts_conn[1], ts_conn[0]
    if v1 not in G: G[v1] = {}
    try:
        ts_id = G[v1][v2]
    except KeyError:
        ts_id = -1 # indicates null value (TS connecting v1 and v2 does not yet exist)
    if ts_id == -1:
        G[v1][v2] = i
    else: # update data structure and file to keep only the lower energy of the two TSs
        if ts_energies[i] < ts_energies[ts_id]:
            vd, vk = ts_id, i # nodes to be deleted/kept, respectively
        else:
            vd, vk = i, ts_id
        G[v1][v2] = vk # update data structure
        print "deleting TS %i, duplicate of TS %i, connect minima %i and %i" % (vd+1, vk+1, v1, v2)
        # subsequent "raw" TS indices need to be consistent with new line numbering
        for j, ts_idx in enumerate(raw_ts_ids[vd+1:]):
            if raw_ts_ids[vd+j+1] > raw_ts_ids[vd]: raw_ts_ids[vd+j+1] -= 1
        # delete line (in-place) corresponding to higher-energy of the two TSs from the ts.data file
        subprocess.call(["sed","-i",str(raw_ts_ids[vd])+"d",tsd_fname])
