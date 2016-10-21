
if (param eq "Te") {
   print "Electron Temperature"
   xcol=1
   ycol=10
   numhead=3
   set style data lines
   leglab="electrons"
   xlab="Position (m)"
   ylab="Temperature (K)" 
} else { if (param eq "ne") {
   print "Electron Density"
   xcol=1
   ycol=3
   numhead=3
   set style data lines
   leglab="electrons"
   xlab="Position (m)"
   ylab="Number Density (m^-3)" 
} else { if (param eq "B") {
   print "Magnetic Field"
   xcol = 1
   ycol = 16
   numhead=3
   set style data lines
   leglab="B_z"
   xlab = "Position (m)"
   ylab = "Magnetic Field (T)"
} else { if (param eq "ni") {
   print "Ion Density"
   xcol=1
   ycol=3
   numhead=3
   set style data lines
   leglab="argon ions"
   xlab="Position (m)"
   ylab="Number Density (m^-3)" 
} else { if (param eq "vi") {
   print "Ion Velocity"
   xcol=1
   ycol=4
   numhead=3
   set style data lines
   leglab="argon v_x"
   xlab="Position (m)"
   ylab="Velocity (m/s)" 
} else { if (param eq "phi") {
   print "Electric Potential"
   xcol = 1
   ycol = 15
   numhead=3
   set style data lines
   leglab="phi"
   xlab = "Position (m)"
   ylab = "Electric Potential (V)"
} else { if (param eq "eps") {
   print "Electron Phase Space"
   xcol = 1
   ycol = 3
   numhead=3
   set style data dots
   leglab="electrons v_x"
   xlab = "Position (m)"
   ylab = "Velocity (m/s)"
} else { if (param eq "ips") {
   print "Ion Phase Space"
   xcol = 1
   ycol = 3
   numhead=3
   set style data dots
   leglab="ions v_x"
   xlab = "Position (m)"
   ylab = "Velocity (m/s)"
} else { if (param eq "evdf") {
   print "Electron Velocity Distribution Function"
   xcol = 3
   ycol = 4
   numhead=2
   set style data lines
   leglab="v_x"
   xlab = "Velocity (m/s)"
   ylab = "evdf"
} else { if (param eq "ivdf") {
   print "Ion Velocity Distribution Function"
   xcol = 3
   ycol = 4
   numhead=2
   set style data lines
   leglab="phi"
   xlab = "Velocity (m/s)"
   ylab = "ivdf"
} else { if (param eq "eedf") {
   print "Electron Energy Distribution Function"
   xcol = 1
   ycol = 2
   numhead=2
   set style data lines
   leglab=""
   xlab = "Velocity (m/s)"
   ylab = "eedf"
} else { if (param eq "iedf") {
   print "Ion Energy Distribution Function"
   xcol = 1
   ycol = 2
   numhead=2
   set style data lines
   leglab=""
   xlab = "Velocity (m/s)"
   ylab = "iedf"
} else {
   print "Unrecognized parameter ".param
   exit
}}}}}}}}}}}}

set xlabel xlab
set ylabel ylab
