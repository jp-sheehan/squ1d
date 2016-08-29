import sys
import math

rdfile = ""
rd = sys.argv[1]
countstr = sys.argv[2]

rdfile = rdfile + rd;

print "Filename: ", rdfile

r = open(rdfile,'r')

tot = 0
count_num = int(countstr) 
errx = 0.0
errx_norm = 0.0
offset = 2
count = 0
m = 9.11e-31
k = 1.38e-23
pi = 3.14158
T = 11600
fxmax = 0.0
fxsum = 0.0


for line in r:
        tot = tot +1
#        if tot>(offset) and count<count_num:
	if tot>(offset):
		vals = line.split()
        	vx = float(vals[2])
                fx = float(vals[3])

                if vx>=0 : 
			fxc = vx*math.sqrt(2*pi*m/(1.38e-23*T))*math.sqrt(m/(2.0*pi*k*T))*math.exp(-m*vx*vx/(2.0*k*T))
                	errx_norm = errx_norm + ((fxc-fx)/fxc)**2
                else:
                	fxc = 0.0

              	errx = errx + (fxc-fx)**2

		fxsum = fxsum+fx

		if fx>fxmax:
			fxmax = fx

		count = count + 1


print "Error: ", math.sqrt(errx/count)
print "Normalized Sum Error: ",  math.sqrt(errx_norm/count)
print "Normalized Max Error: ", math.sqrt(errx/count)/fxmax*100
print "Normalized Mean Error: ", count*math.sqrt(errx/count)/fxsum*100
