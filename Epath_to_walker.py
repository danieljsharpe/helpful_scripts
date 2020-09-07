'''
Read Epath.xxxx files in range (note - list of nodes must be in correct order i.e. corresponding to a real dynamical trajectory, may require reversal of the
original Epath file) and a file "meanwaitingtimes.dat" specifying the mean waiting times for nodes, and write walker.0.x.dat files.
Thus this script produces trajectory files corresponding to the "average" case for k shortest paths
'''

### CHOOSE PARAMS
npaths=1


# read mean waiting times
tau_vals=[]
with open("meanwaitingtimes.dat","r") as tau_f:
    for line in tau_f.readlines():
        tau_vals.append(float(line.split()[0]))

for i in range(npaths):
    node_ids=[] # ordered list of node IDs
    cum_traj_t=[] # accumulated trajectory time values
    cum_t=0.
    curr_Epath_f=open("Epath."+str(i+1).zfill(4),"r")
    for j, line in enumerate(curr_Epath_f.readlines()):
        if j%2!=0: continue # odd-numbered lines correspond to transition states
        line=line.split()
        node_id=int(line[2])
        node_ids.append(node_id)
        cum_traj_t.append(cum_t)
        cum_t+=tau_vals[node_id-1]
    curr_Epath_f.close()
    curr_walker_f=open("walker.0."+str(i)+".dat","w")
    step=0
    for node_id, t in zip(node_ids,cum_traj_t):
        curr_walker_f.write(" %6i    %1i       %5i      %30.10f\n" % (node_id,0,step,t))
        step+=1
    curr_walker_f.close()
