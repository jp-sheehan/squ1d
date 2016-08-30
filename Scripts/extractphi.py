import sys
import math

rdfile = "electronOutput_cField"
wname = "phiout.dat"
exten = ".dat"
start = sys.argv[1]
end = sys.argv[2]

w = open(wname,'w')


for i in range(int(start),int(end)):
	fnum = str(i)
	rname = rdfile+fnum+exten

	print "Filename: ", rname, fnum

	r = open(rname,'r')

	for line in r:
		vals = line.split()

	phi = float(vals[4])
        w.write(fnum+'\t'+str(phi)+'\n')

	print "phi:  ", phi


