''' Print single-column file of n uniformly-distributed random numbers in range x1 to x2
    to file "rand_out.dat" '''

import numpy as np
import sys

n = int(sys.argv[1])
x1 = float(sys.argv[2])
x2 = float(sys.argv[3])

with open("rand_out.dat","w") as outf:
    for i in range(n):
        outf.write("%.6f\n" % np.random.uniform(x1,x2))
