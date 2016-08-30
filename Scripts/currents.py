import sys
import math

rdfile = ""
rd = sys.argv[1]

fnum = sys.argv[2]
start = sys.argv[3]
finish = sys.argv[4]

rdfile = rdfile + rd +fnum+".dat";

print "Filename: ", rdfile

r = open(rdfile,'r')

n = 0
val_avg = 0.0
val1 = 0.0
val2 = 0.0
count = 0

for line in r:
        if n>=int(start) and n<=int(finish):
		vals = line.split()
        	val1 = float(vals[2])
                val2 = float(vals[5])
                val_avg = val1*val2 + val_avg
                count =  count +1
	n = n+1

print "value: ",val_avg/count*1.6e-19
