import sys
import math

rdfile = "ion1Output_pField"
wname = "tempout.dat"
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
			temp = float(vals[9])
			temp1 = float(vals[6])
			temp2 = float(vals[7])
			temp3 = float(vals[8])
			temp4 = float(vals[3])
			temp5 = float(vals[4])
			temp6 = float(vals[5])

        w.write(fnum+'\t'+str(temp)+'\t'+str(temp1)+'\t'+str(temp2)+'\t'+str(temp3)+'\t'+str(temp4)+'\t'+str(temp5)+'\t'+str(temp6)+'\n')
	print "T:  ", temp



