import numpy as np
import sys

fname = sys.argv[1] # single-column file containing data and nothing else
lim = float(sys.argv[2]) # numbers are rounded up to this value if they are below this limit
nlines = int(sys.argv[3])

vec = np.zeros(nlines,dtype=float)

with open(fname,"r") as foo:
    for i, line in enumerate(foo.readlines()):
        x = float(line.split()[0])
        if x<lim: x=lim # round negative numbers to zero
        vec[i] = x

with open("newvals.dat","w") as foo:
    for x in vec:
        foo.write("{:.6e}".format(x)+"\n")
