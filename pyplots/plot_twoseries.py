'''
Simple Python script to plot two data series with two different y-axes
DJS
'''

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.rcParams["text.usetex"]=True # tex font
import matplotlib.pyplot as plt

''' class for plotting two data series with two different y-axes '''
class PlotDouble(object):

    def __init__(self,npts,logvals):
        self.npts = npts
        self.logvals = logvals
        self.xdata = np.zeros(npts,dtype=float)
        self.y1data = np.zeros(npts,dtype=float)
        self.y2data = np.zeros(npts,dtype=float)

    def read_data(self,fname,col1,col2):
        n=0 # count number of lines read
        with open(fname,"r") as foo:
            for line in foo.readlines():
                line = line.split()
                if line[0]=="#": continue # indicates comment line
                self.xdata[n] = float(line[0]) # xdata assumed to be first col
                self.y1data[n] = float(line[col1-1])
                self.y2data[n] = float(line[col2-1])
                n += 1
                if n>=self.npts: break
        if self.logvals:
            self.y1data = np.log10(self.y1data)
            self.y2data = np.log10(self.y2data)

    def plot(self,xlims,ylims1,ylims2,nxticks,nyticks1,nyticks2,figfmt="pdf"):
        # first y axis
        fig, ax1 = plt.subplots(figsize=(10.,7.)) # fig size in inches
        line1 = ax1.plot(self.xdata,self.y1data,marker="o",markersize=10.,color="deeppink",label="$\mathcal{T}_{\mathcal{A}\mathcal{B}}$")
        ax1.set_xlabel("$\mathrm{Inverse\ temperature}\ 1/T$",fontsize=42)
        ax1.set_ylabel("$\log_{10}(\mathrm{MFPT})$",fontsize=42)
        xtick_intvl = float(xlims[1]-xlims[0])/float(nxticks)
        xtick_vals = [xlims[0]+(i*xtick_intvl) for i in range(nxticks+1)]
        ytick1_intvl = float(ylims1[1]-ylims1[0])/float(nyticks1)
        ytick1_vals = [ylims1[0]+(i*ytick1_intvl) for i in range(nyticks1+1)]
        ax1.set_ylim(ylims1)
        ax1.set_xticks(xtick_vals)
        ax1.set_yticks(ytick1_vals)
        xtick_labels = ["$"+"{:.1f}".format(xtick_val)+"$" for xtick_val in xtick_vals]
        ytick1_labels = ["$"+"{:.0f}".format(ytick_val)+"$" for ytick_val in ytick1_vals]
        ax1.set_xticklabels(xtick_labels)
        ax1.set_yticklabels(ytick1_labels)
        # second y axis
        ax2 = ax1.twinx() # second y axis shares same x axis as first y axis
        ax1.tick_params(which="both",direction="out",top=True,labelsize=24)
        line2 = ax2.plot(self.xdata,self.y2data,marker="o",markersize=10.,color="cornflowerblue",label="$\mathcal{J}_{\mathcal{A}\mathcal{B}}$")
        ax1.set_xlim([xlims[0]-(float(xtick_intvl)/10.),xlims[1]+(float(xtick_intvl)/10.)]) # note slight offset on both left and right
        ax2.set_ylabel("$\log_{10}(\mathrm{Reactive\ flux})$",fontsize=42)
        ax2.set_ylim(ylims2)
        ax2.tick_params(direction="out",labelsize=24)
        ytick2_intvl = float(ylims2[1]-ylims2[0])/float(nyticks2)
        ytick2_vals = [ylims2[0]+(i*ytick2_intvl) for i in range(nyticks2+1)]
        ax2.set_yticks(ytick2_vals)
        ytick2_labels = ["$"+"{:.0f}".format(ytick_val)+"$" for ytick_val in ytick2_vals]
        ax2.set_yticklabels(ytick2_labels)
        # legend
        lines = line1+line2
        line_labels = [line.get_label() for line in lines]
        ax1.legend(lines,line_labels,loc="center right",fontsize=28)
        # plot figure and save
        fig.tight_layout()
        plt.savefig("plot_dbl."+figfmt,format=figfmt,bbox_inches="tight")
        plt.show()

if __name__=="__main__":

    ### SET PARAMS ###
    fname = "tot_fluxes.dat" # name of file containing data
    col1 = 3 # index (from 1) of column in file corresponding to first data series
    col2 = 2 # index (from 1) of column in file corresponding to second data series
    npts = 6 # no. of data points
    logvals = True
    ### PLOT PARAMS ###
    xlims = [0.,2.]
    ylims1 = [0.,25.] # for first y axis
    ylims2 = [-25.,0.] # for second y axis
    nxticks = 4
    nyticks1 = 5
    nyticks2 = 5
    logvals = True

    ### RUN ###
    plotobj = PlotDouble(npts,logvals)
    plotobj.read_data(fname,col1,col2)
    plotobj.plot(xlims,ylims1,ylims2,nxticks,nyticks1,nyticks2)
