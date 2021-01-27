# Python script to re-print an Epath file backwards

import sys

Epath_fname = sys.argv[1]

energies = []
sp_idcs = []
n_sp = 0

with open(Epath_fname,"r") as Epath_f:
    for line in Epath_f:
        line = line.split()
        energies.append(float(line[1]))
        sp_idcs.append(int(line[2]))
        n_sp += 1

with open(Epath_fname+".new","w") as out_f:
    for i in range(n_sp):
        out_f.write("  %i   %5.8f    %i\n" % (i+1, energies[-(i+1)], sp_idcs[-(i+1)]))
