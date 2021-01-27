'''
Python script to write min.data.fastest and ts.data.fastest files given original min.data
and ts.data files, and two files mdf_indices_f and tdf_indices_f containing a list of the
minima and transition state indices to appear in the .fastest files.
Can optionally parse the order parameter file also.

Daniel J. Sharpe
Mar 2019
'''

import sys

# read minima and transition state indices, according to the system of the min.data and
# ts.data files, that are to be included in the new .fastest files
def read_sp_idcs(sp_indices_f):
    sp_idcs = []
    with open(sp_indices_f,"r") as spi_f:
        for line in spi_f.readlines():
            sp_idcs.append(int(line))
    return sorted(sp_idcs)

# get new indexing system for minima or transition states
def get_new_sp_idcs_dict(sp_idcs):
    n_sp_idcs_used = 0
    sp_idcs_dict = {}
    for sp_idx in sp_idcs:
        n_sp_idcs_used += 1
        sp_idcs_dict[sp_idx] = n_sp_idcs_used
    return sp_idcs_dict

# write the min.data.fastest.new file
def write_mdf_file(min_idcs, mdf_f, min_idcs_dict):
    mdfnew_f = open("min.data.fastest.new","w")
    with open(mdf_f,"r") as md_f:
        for line_no, line in enumerate(md_f.readlines()):
            if line_no+1 in min_idcs:
                line = line.split()
                mdfnew_f.write("%7.10f     %7.10f     %1i   %6.10f %6.10f %6.10f\n" %
                              (float(line[0]), float(line[1]), int(line[2]), \
                               float(line[3]), float(line[4]), float(line[5])))
    mdfnew_f.close()

# write the ts.data.fastest.new file
def write_tdf_file(ts_idcs, tdf_f, min_idcs_dict):
    tdfnew_f = open("ts.data.fastest.new","w")
    with open(tdf_f,"r") as td_f:
        for line_no, line in enumerate(td_f.readlines()):
            if line_no+1 in ts_idcs:
                line = line.split()
                tdfnew_f.write("%7.10f     %7.10f      %1i  %6i  %6i   %6.10f %6.10f %6.10f\n" %
                              (float(line[0]), float(line[1]), int(line[2]),
                               min_idcs_dict[int(line[3])], min_idcs_dict[int(line[4])],
                               float(line[5]), float(line[6]), float(line[7])))
    tdfnew_f.close()

# write min.A.new and min.B.new files
def write_ab_files(df_app, min_idcs_dict):
    for set_idx in ["A","B"]:
        with open("min."+set_idx+df_app,"r") as min_set_f:
            dummy = min_set_f.readline()
            min_set_idx = int(min_set_f.readline())
        with open("min."+set_idx+".new","w") as min_set_new_f:
            min_set_new_f.write("1\n")
            min_set_new_f.write(str(min_idcs_dict[min_set_idx]))

# write the new order parameter file from the old one
def write_op_file(orderparam_f,min_idcs):
    op_f = open("orderparam_new.dat","w")
    with open(orderparam_f,"r") as old_op_f:
        for line_no, line in enumerate(old_op_f.readlines()):
            if line_no+1 in min_idcs:
                op_f.write(line)
    op_f.close()

if __name__=="__main__":
    mdf_indices_f = sys.argv[1] # file containing list of min indices to be read from min.data
    tdf_indices_f = sys.argv[2] # file containing list of ts indices to be read from ts.data
    if sys.argv[3] != "None": # appendix (e.g. ".all" for original min/ts.data and min.A/B files)
        df_app = sys.argv[3]
    else: # defaults
        df_app = ""
    if len(sys.argv) > 4: orderparam_f = sys.argv[4] # optional
    min_idcs = read_sp_idcs(mdf_indices_f)
    ts_idcs = read_sp_idcs(tdf_indices_f)
    min_idcs_dict = get_new_sp_idcs_dict(min_idcs)
    ts_idcs_dict = get_new_sp_idcs_dict(ts_idcs)
    write_mdf_file(min_idcs,"min.data"+df_app,min_idcs_dict)
    write_tdf_file(ts_idcs,"ts.data"+df_app,min_idcs_dict)
    write_ab_files(df_app,min_idcs_dict)
    if len(sys.argv) > 4: write_op_file(orderparam_f,min_idcs)
