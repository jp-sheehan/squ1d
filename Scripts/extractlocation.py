import sys
import math

pname = "elec1Particles"
exten = ".dat"
filenum = sys.argv[1];
xloc = sys.argv[2];
pfullname = pname + filenum + exten;
writename = "dataextract" + filenum + "_"+ xloc + exten;

fp = open(pfullname,'r');
fw = open(writename,'w');

print "Filename: ", pfullname

print "Number: ", sys.argv[1]

x = float(xloc)
dx = 0.01

x1 = []
x2 = []
v1 = []
v2 = []
v3 = []
en = []

x1n = []
x2n = []
v1n = []
v2n = []
v3n = []
enn = []


next(fp)
next(fp)
next(fp)

for line in fp:
        vals = line.split()
        x1.append(float(vals[0]))
        x2.append(float(vals[1]))
        v1.append(float(vals[2]))
        v2.append(float(vals[3]))
        v3.append(float(vals[4]))
        en.append(float(vals[5]))

for i in range(len(x1)):
	if x1[i]>=(x-0.5*dx) and x1[i]<=(x+0.5*dx):
		x1n.append(x1[i])
		x2n.append(x2[i])
		v1n.append(v1[i])
		v2n.append(v2[i])
		v3n.append(v3[i])
		enn.append(en[i])

for j in range(len(x1n)):
	fw.write(str(x1n[j])+'\t'+str(x2n[j])+'\t'+str(v1n[j]) +'\t'+str(v2n[j])+'\t'+str(v3n[j]) + '\t' + str(enn[j]) + '\n')



