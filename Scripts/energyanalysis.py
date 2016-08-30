import sys
import math

pname1 = "electronParticles"
pname2 = "ionParticles"
fname1 = "electronOutput_pField"
exten = ".dat"
filenum = sys.argv[1];
pfullname1 = pname1 + filenum + exten;
pfullname2 = pname2 + filenum + exten;
ffullname1 = fname1 + filenum + exten;

fp1 = open(pfullname1,'r');
fp2 = open(pfullname2,'r');
ff1 = open(ffullname1,'r');

print "Filename 1: ", pfullname1
print "Filename 2: ", pfullname2
print "Filename 3: ", ffullname1

print "Number: ", sys.argv[1]

en_average = 0.0
error = 0.0
q =  0.049087
m = 0.049087
dx = 2*3.14159/32;

fw = open('energy.dat','w')

enke_1 = []
enke_2 = []
enE = []

avg_enke = 0.0
avg_entot = 0.0
avg_enE = 0.0
count= 0

for line in fp1:
        vals = line.split()
        enke_1.append(float(vals[5]))

for line in fp2:
        vals = line.split()
        enke_2.append(float(vals[5]))

for line in ff1:
        vals = line.split()
        enE.append(float(vals[2]))

for i in range(len(enke_1)):
	avg_enke = avg_enke + enke_1[i] + enke_2[i];

for i in range(len(enE)):
	avg_enE = avg_enE + (0.5)*enE[i]*enE[i]*dx;


print "Field Energy:\t",avg_enE 
print "Kinetic Energy:\t",avg_enke
print "Total Energy:\t",avg_enke+avg_enE
