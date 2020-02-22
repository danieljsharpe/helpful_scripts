'''
Read min.data and ts.data files and write min.data.new and ts.data.new files with no frequencies (ie frequencies replaced by dummy values 0/1)
'''

newmd_f = open("min.data.new","w")
with open("min.data","r") as oldmd_f:
    for line in oldmd_f.readlines():
        line=line.split()
        newmd_f.write("%.10f   %1.10f    %1i  %.10f  %.10f  %.10f\n" % \
                      (float(line[0]),0.,1,float(line[3]),float(line[4]),float(line[5])))
newmd_f.close()

newtsd_f = open("ts.data.new","w")
with open("ts.data","r") as oldtsd_f:
    for line in oldtsd_f.readlines():
        line=line.split()
        newtsd_f.write("%.10f   %1.10f   %1i   %6i  %6i   %.10f  %.10f  %.10f\n" % \
                       (float(line[0]),0.,1,int(line[3]),int(line[4]),float(line[5]),float(line[6]),float(line[7])))
newtsd_f.close()
