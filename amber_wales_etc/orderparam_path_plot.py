'''
Python script to plot evolution of 2 order parameters along an energy surface
Input arguments:
1. Epath file start index  2. Epath file end index  3. order parameter file #1  4. order parameter file #2
5. max path length (of minima only)
6. / 7. min / max values of order param 1 (for plotting)
8. / 9. min / max values of order param 2 (for plotting)
10. mode. If =1: plot 2 order params vs energy. If=2: plot 3 order params and colour by energy
11. (only if mode=2) order parameter file #3
12. / 13. (only if mode=2) min / max values of order param 3 (for plotting)
Need to have Epath files Epath.x - Epath.y in your current working directory
order parameter files must be single-column

Daniel J. Sharpe
Mar 2019
'''

import sys
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def read_op_vals(op_fname):
    op_vals = []
    with open(op_fname,"r") as op_f:
        for line in op_f.readlines():
            op_val = float(line)
            if op_val==0.: op_val = 1.0E-8 # so that trim_zeros() doesn't mess up the array sizes
            op_vals.append(op_val)
    return op_vals

def get_path_info(path_no_first,path_no_last,max_path_len,op1_vals,op2_vals,op3_vals,plot_mode):
    n_paths = path_no_last-path_no_first+1
    path_info_mtx = np.zeros((n_paths,3,max_path_len),dtype=float)
    for path_idx, path_no in enumerate(range(path_no_first,path_no_last+1)):
        with open("Epath."+str(path_no),"r") as ep_f:
            min_no = 0
            for step_no, line in enumerate(ep_f.readlines()):
                if (step_no+1)%2==0: continue # ignore TSs
                line = line.split()
                min_id = int(line[2])
                if plot_mode==1: path_info_mtx[path_idx,0,min_no] = float(line[1]) # energy value
                elif plot_mode==2: path_info_mtx[path_idx,0,min_no] = op3_vals[min_id-1]
                path_info_mtx[path_idx,1,min_no] = op1_vals[min_id-1] # order param 1 value
                path_info_mtx[path_idx,2,min_no] = op2_vals[min_id-1] # order param 2 value
                min_no += 1
    return path_info_mtx

def plot_op_evoln(path_info_mtx):
    my_cmap = plt.cm.brg(np.linspace(0,1,np.shape(path_info_mtx)[0]))
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    for i in range(np.shape(path_info_mtx)[0]):
        xdata = np.trim_zeros(path_info_mtx[i,1,:],"b")
        ydata = np.trim_zeros(path_info_mtx[i,2,:],"b")
        zdata = np.trim_zeros(path_info_mtx[i,0,:],"b")
        ax.plot(xdata,ydata,zdata,color=my_cmap[i])
    ax.set_xlabel("alpha helix")
    ax.set_ylabel("random coil")
    ax.set_zlabel("beta sheet")
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("w")
    ax.yaxis.pane.set_edgecolor("w")
    ax.zaxis.pane.set_edgecolor("w")
#    plt.axis("off")
#    plt.grid(b=None)
    ax.grid(False)
    plt.show()

if __name__=="__main__":
    path_no_first = int(sys.argv[1])
    path_no_last = int(sys.argv[2])
    op1_fname = sys.argv[3]
    op2_fname = sys.argv[4]
    max_path_len = int(sys.argv[5])
    op1_min, op1_max = float(sys.argv[6]), float(sys.argv[7])
    op2_min, op2_max = float(sys.argv[8]), float(sys.argv[9])
    plot_mode = int(sys.argv[10])
    if plot_mode==2:
        op3_fname = sys.argv[11]
        op3_min, op3_max = float(sys.argv[12]), float(sys.argv[13])

    op1_vals = read_op_vals(op1_fname)
    op2_vals = read_op_vals(op2_fname)
    if plot_mode==2: op3_vals = read_op_vals(op3_fname)
    else: op3_vals = None

    path_info_mtx = get_path_info(path_no_first,path_no_last,max_path_len,op1_vals,op2_vals,op3_vals,plot_mode)

    plot_op_evoln(path_info_mtx)
