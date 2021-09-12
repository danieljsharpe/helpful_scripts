import numpy as np
import sys

fname = sys.argv[1] # single-column file containing data and nothing else
nlines = int(sys.argv[2])

vec = np.zeros(nlines,dtype=float)

with open(fname,"r") as foo:
    for i, line in enumerate(foo.readlines()):
        x = np.log10(float(line.split()[0])) # log ten of data
        vec[i] = x

with open("newvals.dat","w") as foo:
    for x in vec:
        foo.write("{:.6e}".format(x)+"\n")
