import sys
import math

rdfile = "elec1Particles"
wname = "timehistory.dat"
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

	x = float(vals[0])
	y = float(vals[1])
	v1 = float(vals[2])
	v2 = float(vals[3])
	v3 = float(vals[4])
	en = float(vals[5])
        w.write(str(x)+'\t'+str(v1)+'\t'+str(v2)+'\t'+str(v3)+'\t'+str(en)+'\n')

	print "en:  ", en


