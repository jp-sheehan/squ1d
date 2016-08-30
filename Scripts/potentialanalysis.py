import sys
import math

name = "PICfield"
exten = ".dat"
filenum = sys.argv[1];
fullname = name + filenum + exten;
fr = open(fullname,'r');

print "Filename: ", sys.argv[0]
print "Number: ", sys.argv[1]

average = 0.0
error = 0.0
count = 0

for line in fr:
        vals = line.split()
	phi = float(vals[6])
        x = float(vals[0])
        phiexact = -17920.7*x*x/2.0 + 17920.7/2.0*x 

        average = average + phi
        error = error + math.fabs(phi-phiexact)
        count = count + 1

print "Error:   ", error/count
print "Average:   ", average/count




