import sys
import math

name1 = "electronParticles"
name2 = "ionParticles"
exten = ".dat"
filenum = sys.argv[1];
fullname1 = name1 + filenum + exten;
fullname2 = name2 + filenum + exten;
fullnamein1 = name1 + "0" + exten;
fullnamein2 = name2 + "0" + exten;

fr1 = open(fullname1,'r');
fr2 = open(fullname2,'r');
fin1 = open(fullnamein1,'r');
fin2 = open(fullnamein2,'r');

print "Filename 1: ", fullname1
print "Filename 2: ", fullname2
print "Number: ", sys.argv[1]

en_average = 0.0
error = 0.0
q =  0.049087
m = 0.049087

fw = open('energy.dat','w')

enphi_in_1 = []
enke_in_1 = []
entot_in_1 = []

enphi_in_2 = []
enke_in_2 = []
entot_in_2 = []

enphi_1 = []
enke_1 = []
entot_1 = []

enphi_2 = []
enke_2 = []
entot_2 = []

avg_enphi_in = 0.0
avg_entot_in = 0.0
avg_enke_in = 0.0
count= 0

for line in fin1:
        vals = line.split()
        enphi_in_1.append(float(vals[6]))
        enke_in_1.append(float(vals[5]))
        entot_in_1.append(float(vals[7]))

for line in fin2:
        vals = line.split()
        enphi_in_2.append(float(vals[6]))
        enke_in_2.append(float(vals[5]))
        entot_in_2.append(float(vals[7]))


avg_enphi = 0.0
avg_entot = 0.0
avg_enke = 0.0
count= 0

for line in fr1:
        vals = line.split()
        enphi_1.append(float(vals[6]))
        enke_1.append(float(vals[5]))
        entot_1.append(float(vals[7]))

for line in fr2:
        vals = line.split()
        enphi_2.append(float(vals[6]))
        enke_2.append(float(vals[5]))
        entot_2.append(float(vals[7]))


for i in range(len(enphi_1)):
	avg_enphi_in = avg_enphi_in + abs(enphi_in_1[i]) + abs(enphi_in_2[i]);
	avg_enke_in = avg_enke_in + enke_in_1[i] + enke_in_2[i];
	avg_entot_in = avg_entot_in + entot_in_1[i] + entot_in_2[i];

	avg_enphi = avg_enphi + abs(enphi_1[i]) + abs(enphi_2[i]);
	avg_enke = avg_enke + enke_1[i] + enke_2[i];
	avg_entot = avg_entot + entot_1[i] + entot_2[i];

print "Field Energy:\t",avg_enphi 
print "Kinetic Energy:\t",avg_enke
print "Total Energy:\t", avg_entot
print "Initial Energy:\t", avg_enphi_in

