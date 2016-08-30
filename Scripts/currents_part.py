import sys
import math

rdfile = ""
rd = sys.argv[1]

fnum = sys.argv[2]
x1 = sys.argv[3]
x2 = sys.argv[4]

rdfile = rdfile + rd +fnum+".dat";

print "Filename: ", rdfile

r = open(rdfile,'r')

n = 0
start = 3
val_sum =0.0
val_avg = 0.0
val1 = 0.0
val2 = 0.0
count = 0

for line in r:
        if n>=int(start):
		vals = line.split()
		if float(vals[0])>=float(x1) and float(vals[0])<=float(x2):
        		val1 = float(vals[4])
                	val_sum = val1 + val_sum
                	count =  count +1
	n = n+1

print "value: ",val_sum*1.6e-19*1e9/(float(x2)-float(x1)), val_sum/count, count
