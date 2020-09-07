'''
parse a file corresponding to node data and if the value is greater than or equal to a specified threshold (default 0.), add the node to a list
These nodes are written to file nodelist.dat
'''

import sys
from math import isnan

valsfile = sys.argv[1]
try:
    threshval = float(sys.argv[2])
except IndexError:
    threshval = 0.

nodelist=[]
with open(valsfile,"r") as vals_f:
    for i, line in enumerate(vals_f.readlines()):
        val=float(line.split()[0])
        if isnan(val): continue
        if val>=threshval: nodelist.append(i+1)
with open("nodelist.dat","w") as nodes_f:
    for node in nodelist:
        nodes_f.write("%i\n" % node)
print "written %i node IDs to file" % len(nodelist)
