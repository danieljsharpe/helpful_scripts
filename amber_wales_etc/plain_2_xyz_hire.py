'''
Python script to parse a plain coordinates file (extractedmin/extractedts -style) of a HiRE-RNA sequence
'''

import numpy as np
import sys

# the base sequence corresponding to the coordinates file being processed
hire_seq = "AGGGTTAGGGTTAGGGTTAGGG"


# no. of beads used to represent each nucleotide
nbeads_dict = {"A": 7, "G": 7, "T": 6, "C": 6}

# NB the plain coords file should have no blank lines
def read_plain_coords_file(fname,natoms):
    coords = np.zeros((natoms,3),dtype=float)
    with open(fname,"r") as coordsf:
        for i, line in enumerate(coordsf.readlines()):
            line = line.split()
            curr_coords = np.array([float(line[0]),float(line[1]),float(line[2])])
            coords[i,:] = curr_coords
    return coords

def write_xyz_file(fname,coords,base_ids):
    with open(fname,"w") as outf:
        outf.write("%4i\n" % len(base_ids))
        outf.write("\n")
        for i in range(len(base_ids)):
            outf.write("%4s   %6.6f  %6.6f  %6.6f\n" % (base_ids[i], coords[i,0], coords[i,1], coords[i,2]))

if __name__=="__main__":
    infile = sys.argv[1]
    outfile = infile+".xyz"
    natoms = 0
    base_ids = []
    for base in hire_seq:
        natoms += nbeads_dict[base]
        if base_ids:
            base_ids.extend([base]*nbeads_dict[base])
        else: # first (5') nucleotide omits phosphate, so has one less bead
            base_ids.extend([base]*(nbeads_dict[base]-1))
    coords = read_plain_coords_file(infile,natoms)
    write_xyz_file(outfile,coords,base_ids)
