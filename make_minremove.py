removemin = []
nmin = 0

with open("disconnected","r") as disf:
    for line in disf.readlines():
        line = line.split()
        removemin.append(str(line[2]))
        nmin += 1

# NB also need a file ts.remove that is just one line with "0"
with open("min.remove","w") as mrf:
    mrf.write(str(nmin)+"\n")
    for minimum in removemin:
        mrf.write(str(minimum)+"\n")
