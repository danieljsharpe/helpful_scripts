'''
parse a min.data file and categorise the minima depending on their energy (useful for single-funnel systems). Can also specify separate category of absorbing nodes
'''

### SET THE INTERVALS
intvl_list = [-427.,-439.]
### SET ABSORBING NODES
abs_nodes = [2,14804,14934,15899,21148,31232,42054,47758,62421,64682,64972,65290,66696,66810,67215,67987,68179]

categories=[]
with open("min.data","r") as md_f:
    for line in md_f.readlines():
        energy=float(line.split()[0])
        idx=0
        while idx<len(intvl_list) and intvl_list[idx]>energy: idx+=1
        categories.append(idx)
if abs_nodes:
    for node_id in abs_nodes:
        categories[node_id-1]=len(intvl_list)+1
with open("superstates.dat","w") as cats_f:
    for cat in categories:
        cats_f.write("%i\n" % cat)
