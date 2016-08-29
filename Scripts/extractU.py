import sys
import math

rdfile = "ion1Output_pField"
wname = "Uout.dat"
exten = ".dat"
start = sys.argv[1]
end = sys.argv[2]

w = open(wname,'w')

i = 0
loc = 5

for i in range(int(start),int(end)):
	fnum = str(i)
	rname = rdfile+fnum+exten

	print "Filename: ", rname, fnum

	r = open(rname,'r')

        i=0  

	for line in r:
                i=i+1
		vals = line.split()
                if(i==loc):
			temp = float(vals[3])

        w.write(fnum+'\t'+str(temp)+'\n')
	print "T:  ", temp



