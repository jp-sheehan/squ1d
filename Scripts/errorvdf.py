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
err = 0.0
err_norm = 0.0
errx = 0.0
errx_norm = 0.0
offset = 2
count = 0
m = 9.11e-31
k = 1.38e-23
pi = 3.14158
T = 11600
ftotmax = 0.0
fxmax = 0.0
ftotsum = 0.0
fxsum = 0.0


for line in r:
        tot = tot +1
#        if tot>(offset) and count<count_num:
	if tot>(offset):
		vals = line.split()
        	vtot = float(vals[0])
                ftot = float(vals[1])
                ftotc = math.sqrt((m/(2.0*pi*k*T))**3)*4.0*pi*vtot*vtot*math.exp(-m*vtot*vtot/(2*k*T))
                err = err + (ftotc-ftot)**2
                err_norm = err_norm + ((ftotc-ftot)/ftotc)**2

        	vx = float(vals[2])
                fx = float(vals[3])
                fxc = math.sqrt(m/(2.0*pi*k*T))*math.exp(-m*vx*vx/(2.0*k*T))
                errx = errx + (fxc-fx)**2
                errx_norm = errx_norm + ((fxc-fx)/fxc)**2

		ftotsum = ftotsum+ftot
		fxsum = fxsum+fx

		if ftot>ftotmax:
			ftotmax = ftot
		if fx>fxmax:
			fxmax = fx

		count = count + 1


print "Error: ", math.sqrt(err/count),"\t" , math.sqrt(errx/count)
print "Normalized Sum Error: ", math.sqrt(err_norm/count),"\t", math.sqrt(errx_norm/count)
print "Normalized Max Error: ", math.sqrt(err/count)/ftotmax*100,"\t", math.sqrt(errx/count)/fxmax*100
print "Normalized Mean Error: ", count*math.sqrt(err/count)/ftotsum*100,"\t", count*math.sqrt(errx/count)/fxsum*100
