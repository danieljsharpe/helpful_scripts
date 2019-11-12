'''
Python script to parse the file "disconnected" prduced by disconnectionDPS, and rewrite the
min.data.dummy, ts.data.dummy, min.pos.dummy and ts.pos.dummy files
not including those minima that are disconnected

Daniel J. Sharpe
Feb 2019
'''

import numpy as np

def read_sp_files(spfile,ftype,disconn_id=None,bad_ts=None):
    sp_info = []
    with open(spfile,"r") as spf:
        for sp_id, line in enumerate(spf.readlines()):
            line = line.split()
            if ftype==0: # reading min.data.dummy for minima energies
                if sp_id+1 not in disconn_id:
                    sp_info.append(float(line[0]))
            elif ftype==1: # reading min.pos.dummy for minima positions
                if sp_id+1 not in disconn_id:
                    sp_info.append([float(x) for x in line])
            elif ftype==2: # reading ts.pos.dummy for ts positions
                if sp_id+1 not in bad_ts:
                    sp_info.append([float(x) for x in line])
    return sp_info

def get_new_minid(disconn_id,min_id):
    return disconn_id.index(min_id)

def read_tsdata_file(tsdfile,disconn_id):
    ts_info = []
    bad_ts = []
    with open("ts.data.dummy","r") as tsdf:
        for ts_id, line in enumerate(tsdf.readlines()):
            line = line.split()
            conn1, conn2 = int(line[3]), int(line[4])
            diff_bad = [0,0]
            ts_is_bad = False
            for bad_min in disconn_id:
                if bad_min < conn1: diff_bad[0] += 1
                if bad_min < conn2: diff_bad[1] += 1
                if (bad_min == conn1) or (bad_min == conn2):
                    bad_ts.append(ts_id+1)
                    ts_is_bad = True
                    break
            if not ts_is_bad:
                ts_info.append([float(line[0]),conn1-diff_bad[0],conn2-diff_bad[1]])
    return ts_info, bad_ts

def write_data_files(datafile,sp_data,ftype):
    with open(datafile,"w") as df:
        if ftype==0: # write min.data -style file
            for entry in sp_data:
                df.write("%f %f %i %f %f %f\n" % (entry,1.,1,1.,1.,1.))
        elif ftype==1: #write ts.data -style file
            for entry in sp_data:
                df.write("%f %f %i %i %i %f %f %f\n" %
                         (entry[0],0.,1,entry[1],entry[2],1.,1.,1.))

def write_pos_files(posfile,posns):
    with open(posfile,"w") as posf:
        for pos in posns:
            for x in pos:
                posf.write("  %f" % x)
            posf.write("\n")

if __name__ == "__main__":
    # read minima that are disconnected
    disconn_id = []
    with open("disconnected","r") as disconnf:
        for line in disconnf.readlines():
            line = line.split()
            disconn_id.append(int(line[2]))
    # read in data
    min_data = read_sp_files("min.data.dummy",0,disconn_id)
    min_pos = read_sp_files("min.pos.dummy",1,disconn_id)
    ts_data, bad_ts = read_tsdata_file("ts.data.dummy",disconn_id)
    ts_pos = read_sp_files("ts.pos.dummy",2,bad_ts=bad_ts)
    # write out data
    write_data_files("min.data.dummy.removed",min_data,0)
    write_data_files("ts.data.dummy.removed",ts_data,1)
    write_pos_files("min.pos.dummy.removed",min_pos)
    write_pos_files("ts.pos.dummy.removed",ts_pos)
