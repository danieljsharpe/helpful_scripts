'''
Script to process a pdb file for a DNA/RNA structure and modify to give template RNA/DNA, respectively.
Then load into Amber Leap to get proper pdb file.
For use with OL15 FF.

Usage:
python DNA_2_RNA.py foo-in.pdb bar-out.pdb

Daniel J. Sharpe
September 2018
'''

import sys


pdb_fin = sys.argv[1]
pdb_fout = sys.argv[2]

DA_2_RA, RA_2_DA = False, False
with open(pdb_fin,"r") as fin:
    fout = open(pdb_fout,"w")
    for line_no, line in enumerate(fin.readlines()):
        try:
            line = line.split()
            line[2]
        except IndexError: # line corresponds to remark in pdb file
            continue
        if line_no == 1: # use first line to determine if we have DNA or RNA
            if line[3][0] == "D":
                DA_2_RA = True
            elif line[3][0] in ["G","C","A","U"]:
                RA_2_DA = True
            if not DA_2_RA and not RA_2_DA:
                quit("Error: appear to have neither DNA nor RNA")
        if DA_2_RA and line[2] == "H2''":
            line[2] = "O2'"
        elif RA_2_DA and line[2] == "O2'":
            line[2] = "H2''"
        elif RA_2_DA and line[2] == "HO2'":
            continue
        if DA_2_RA and line[3][1] == "T":
            res = list(line[3])
            res[1] = "U"
            line[3] = "".join(res)
            if line[2][0:2] == "H7" or line[2] == "C7": continue
        elif RA_2_DA and line[3][1] == "U":
            res = list(line[3])
            res[1] = "T"
            line[3] = "".join(res)
            if line[2] == "H5": line[2] = "C7"
        if DA_2_RA: line[3] = line[3][1:]
        elif RA_2_DA: line[3][0] = "D"
        fout.write("%4s   %4s %-4s %3s    %2s     %7s %7s %7s  %.2f  %.2f\n" \
                       % (line[0], line[1], line[2], line[3], line[4], line[5], \
                          line[6], line[7], 1.00, 0.00))
    fout.write("TER\nEND")
    fout.close()
