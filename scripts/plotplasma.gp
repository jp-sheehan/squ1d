
set print "-"
set terminal x11
set output

load "~/squ1d/scripts/setplotparams.gp"

if (param eq "Te") {
   filename=sprintf("%selectronOutput_pField%d.dat",dir,fnum)
} else { if (param eq "ne") {
   filename=sprintf("%selectronOutput_pField%d.dat",dir,fnum)
} else { if (param eq "B") {
   filename=sprintf("%selectronOutput_pField%d.dat",dir,fnum)
} else { if (param eq "ni") {
   filename=sprintf("%s%sOutput_pField%d.dat",dir,ionspec,fnum)
} else { if (param eq "vi") {
   filename=sprintf("%s%sOutput_pField%d.dat",dir,ionspec,fnum)
} else { if (param eq "phi") {
   filename=sprintf("%selectronOutput_pField%d.dat",dir,fnum)
} else { if (param eq "eps") {
   filename=sprintf("%selectronParticles%d.dat",dir,fnum)
} else { if (param eq "ips") {
   filename=sprintf("%s%sParticles%d.dat",dir,ionspec,fnum)
} else { if (param eq "evdf") {
   filename=sprintf("%selectrontotalvdf%d.dat",dir,fnum)
} else { if (param eq "ivdf") {
   filename=sprintf("%s%stotalvdf%d.dat",dir,ionspec,fnum)
} else { if (param eq "eedf") {
   filename=sprintf("%selectrontotalvdf%d.dat",dir,fnum)
} else { if (param eq "iedf") {
   filename=sprintf("%s%stotalvdf%d.dat",dir,ionspec,fnum)
} else {
   print "Unrecognized parameter ".param
   exit
}}}}}}}}}}}}
print "from ".filename

plot filename\
   using xcol:ycol\
   every ::numhead\
   title leglab
pause -1 "Hit return to quit"
