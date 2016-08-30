import sys
import math

name = "electronOutput_cField"
exten = ".dat"
filenum = sys.argv[1];
fullname = name + filenum + exten;
fullnamein = name + "0" + exten;

fr = open(fullname,'r');
fin = open(fullnamein,'r');

print "Filename: ", sys.argv[0]
print "Number: ", sys.argv[1]

en_average = 0.0
error = 0.0
q = 0.049087
m = 0.049087

fw = open('energy.dat','w')
phi_in = []

avg_phi_in = 0.0
count= 0
for line in fin:
        vals = line.split()
        rho = float(vals[2])
        phi_in.append(float(vals[7]))
        avg_phi_in = avg_phi_in + 0.5*q*rho*phi_in[count]

        if count==0:
		avg_phi_in = avg_phi_in - 0.5*q*rho*phi_in[count]
        if count==1:
		avg_phi_in = avg_phi_in - 0.5*q*rho*phi_in[count]

        count = count+1


avg_phi = 0.0
avg_en = 0.0
count = 0
for line in fr:
        vals = line.split()
	en = float(vals[6])
        x = float(vals[0])
        phi = float(vals[7])
        rho = float(vals[2])

        avg_en = avg_en + rho*en;
        avg_phi = avg_phi + 0.5*q*rho*phi;

        en_phi = rho*(phi-phi_in[count])*q
        en_vel = rho*en

        if count==0:
		avg_en = avg_en - rho*en
		avg_phi = avg_phi - 0.5*q*rho*phi
        if count==1:
		avg_en = avg_en - rho*en
		avg_phi = avg_phi - 0.5*q*rho*phi

        fw.write(str(x))
        fw.write("\t")
        fw.write(str(en_vel))
        fw.write("\t")
        fw.write(str(en_phi))
        fw.write("\n")
        count = count+1

print "Field Energy:\t",avg_phi 
print "Kinetic Energy:\t",avg_en
print "Total Energy:\t", avg_en+avg_phi
print "Initial Energy:\t", avg_phi_in

