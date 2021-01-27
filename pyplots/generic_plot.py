'''
A general and adaptable Python script to plot line graphs
DJS
'''

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.rcParams["text.usetex"]=True # tex font
import matplotlib.pyplot as plt

class PlotMulti(object):

    def __init__(self,nseries,nent):
        self.nseries=nseries
        self.nent=nent
        self.data = np.zeros((nseries,nent),dtype=float) # array in which data is stored
        self.xvals = np.zeros(nent,dtype=int) # x values are read from one file and are assumed to be the same for all data series

    def read_data(self,froot,snames,col):
        for k, sname in enumerate(snames):
            foo=open(froot+"."+sname+".dat","r")
            for i in range(self.nent):
                line = foo.readline()
                self.data[k,i] = float(line.split()[col-1])
                if k==0: self.xvals[i] = float(line.split()[0])
            foo.close()

    def plot(self,xlims,ylims,slabels,nxticks,nyticks):
        colors=plt.cm.Blues(np.linspace(0.4,1,self.nseries)) # colormap. Note offset to avoid more transparent region of colormap
        plt.figure(figsize=(10.,7.)) # fig size in inches
        for k in range(self.nseries):
            plt.plot(self.xvals,self.data[k],linewidth=6,label="$"+"{:.2f}".format(slabels[k])+"$",color=colors[k])
        plt.xlabel("$\mathrm{No.\ of\ transition\ flux}$-$\mathrm{paths}$",fontsize=42)
        plt.ylabel("$\mathrm{Cumulative\ relative\ flux}$",fontsize=42)
        ax = plt.gca()
        # axis limits; tick values and labels
        ax.tick_params(direction="out",labelsize=24,top=True,right=True)
        xtick_intvl = float(xlims[1]-xlims[0])/float(nxticks)
        assert(xtick_intvl.is_integer())
        xtick_vals = [xlims[0]+(i*xtick_intvl) for i in range(nxticks+1)]
        ytick_intvl = float(ylims[1]-ylims[0])/float(nyticks)
        ytick_vals = [ylims[0]+(i*ytick_intvl) for i in range(nyticks+1)]
        ax.set_xlim([-1.*float(xtick_intvl)/20.,xlims[1]]) # note slight offset
        ax.set_ylim([ylims[0],ylims[1]+1.*float(ytick_intvl)/20.]) # note slight offset
        ax.set_xticks(xtick_vals)
        ax.set_yticks(ytick_vals)
        xtick_labels = ["$"+"{:d}".format(int(xtick_val))+"$" for xtick_val in xtick_vals]
        ytick_labels = ["$"+"{:.1f}".format(ytick_val)+"$" for ytick_val in ytick_vals]
        ax.set_xticklabels(xtick_labels)
        ax.set_yticklabels(ytick_labels)
        # annotate: put labels on right hand side of figure adjacent to final data points, for each series
        for k in range(self.nseries):
            ax.annotate("$"+"{:.2f}".format(slabels[k])+"$",(0.85,0.95-(0.05*(6-k))), \
                        xycoords=("axes fraction","data"),color=colors[k],fontsize=24)
        # save and show plot
        plt.tight_layout()
        plt.savefig("plot_series.pdf",fmt="pdf",bbox_inches="tight")
        plt.show()

if __name__=="__main__":

    ### PARAMS ###
    nseries = 6 # no. of files to read from (=no. of data series)
    froot = "fpp_properties" # root file name. File names are froot.snames[n].dat for all n
    snames = ["010","050","100","150","175","200"] # series names in files
    nent = 20 # number of entries in each data file
    col = 4 # index (from 1) of column from which to read
    logvals = False
    ### PLOT PARAMS ###
    slabels = [0.1,0.5,1.,1.5,1.75,2.] # series labels for plot
    xlims = [0,15]
    ylims = [0.3,1.]
    nxticks = 3
    nyticks = 8

    ### RUN ###
    assert(nseries==len(snames) and len(snames)==len(slabels))

    plotobj = PlotMulti(nseries,nent)
    plotobj.read_data(froot,snames,col)

    # calc accumulated probabilities for path action
    plotobj.data = np.exp(-1.*plotobj.data) # since: path_prob = exp(-1.*path_cost)
    for k in range(nseries):
        for i in range(1,nent):
            plotobj.data[k,i] += plotobj.data[k,i-1]

    # calc the accumulated *relative* contributions
    for k in range(nseries):
        plotobj.data[k,:] *= 1./plotobj.data[k,-1]

    if logvals: plotobj.data = np.log10(plotobj.data)

    plotobj.plot(xlims,ylims,slabels,nxticks,nyticks)
