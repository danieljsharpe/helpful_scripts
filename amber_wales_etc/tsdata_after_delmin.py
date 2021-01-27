'''
Write ts.data.new file, deleting transition states that connect a minimum with ID greater than the argument
'''

import sys

max_min_id = int(sys.argv[1])

tdnew_f = open("ts.data.new","w")
with open("ts.data","r") as td_f:
    for line in td_f.readlines():
        line = line.split()
        if not ((int(line[3]) > max_min_id) or (int(line[4]) > max_min_id)):
            tdnew_f.write("%7.10f     %7.10f      %1i  %6i  %6i   %6.10f %6.10f %6.10f\n" %
                         (float(line[0]), float(line[1]), int(line[2]),
                          int(line[3]), int(line[4]),
                          float(line[5]), float(line[6]), float(line[7])))
tdnew_f.close()
