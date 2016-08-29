import sys
import math

pname = "elec1Particles"
exten = ".dat"
filenum = sys.argv[1];
loc = sys.argv[2];
pfullname = pname + filenum + exten;
writename = "pdf" + filenum + exten;

fp = open(pfullname,'r');
fw = open(writename,'w');

print "Filename: ", pfullname

print "Number: ", sys.argv[1]

numvels = 333

maxvel  =  0.03
delvel  =  2*maxvel/numvels
x = float(loc)
dx = 0.1

pdf   =  [0]*numvels
pdf1  =  [0]*numvels
pdf2  =  [0]*numvels
pdf3  =  [0]*numvels

v1pdf = []
v2pdf = []
v3pdf = []
vmagpdf = []

for i in range(len(pdf)):
	v1pdf.append(delvel*((i-numvels/2)+0.5))
	v2pdf.append(delvel*((i-numvels/2)+0.5))
	v3pdf.append(delvel*((i-numvels/2)+0.5))
	vmagpdf.append(delvel*((i-numvels/2)+0.5))

x1 = []
v1 = []
v2 = []
v3 = []
vmag = []

next(fp)
next(fp)
next(fp)

for line in fp:
        vals = line.split()
        x1.append(float(vals[0]))
        v1.append(float(vals[2]))
        v2.append(float(vals[3]))
        v3.append(float(vals[4]))

for i in range(len(v1)):
	vmag.append(math.sqrt(v1[i]*v1[i] + v2[i]*v2[i] + v3[i]*v3[i]))

	if x1[i]>=(x-0.5*dx) and x1[i]<=(x+0.5*dx):
		for j in range(len(pdf)):
			if vmag[i]>=(vmagpdf[j]-0.5*delvel) and vmag[i]<=(vmagpdf[j]+0.5*delvel):
				pdf[j] = pdf[j] +1
			if v1[i]>=(v1pdf[j]-0.5*delvel) and v1[i]<=(v1pdf[j]+0.5*delvel):
				pdf1[j] = pdf1[j] +1
			if v2[i]>=(v2pdf[j]-0.5*delvel) and v2[i]<=(v2pdf[j]+0.5*delvel):
				pdf2[j] = pdf2[j] +1
			if v3[i]>=(v3pdf[j]-0.5*delvel) and v3[i]<=(v3pdf[j]+0.5*delvel):
				pdf3[j] = pdf3[j] +1


for j in range(len(pdf)):
	fw.write(str(vmagpdf[j])+'\t'+str(pdf[j]/(len(vmag)*delvel))+'\t'+str(pdf1[j]/(len(vmag)*delvel))+'\t'+str(pdf2[j]/(len(vmag)*delvel))+'\t'+str(pdf3[j]/(len(vmag)*delvel))+'\n')



