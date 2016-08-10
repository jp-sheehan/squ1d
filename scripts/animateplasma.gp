set print "-"
set terminal gif animate delay 100
set output sprintf("%s%s.gif",dir,param)

load "~/squ1d/scripts/setplotparams.gp"

if (param eq "Te") {
   filename_base=sprintf("%selectronOutput_%sField",dir,loc)
} else { if (param eq "ne") {
   filename_base=sprintf("%selectronOutput_%sField",dir,loc)
} else { if (param eq "B") {
   filename_base=sprintf("%selectronOutput_%sField",dir,loc)
} else { if (param eq "ni") {
   filename_base=sprintf("%s%sOutput_%sField",dir,ionspec,loc)
} else { if (param eq "vi") {
   filename_base=sprintf("%s%sOutput_%sField",dir,ionspec,loc)
} else { if (param eq "phi") {
   filename_base=sprintf("%selectronOutput_%sField",dir,loc)
} else { if (param eq "eps") {
   filename_base=sprintf("%selectronParticles",dir)
} else { if (param eq "ips") {
   filename_base_base=sprintf("%s%sParticles",dir,ionspec)
} else { if (param eq "evdf") {
   filename_base=sprintf("%selectrontotalvdf",dir)
} else { if (param eq "ivdf") {
   filename_base=sprintf("%s%stotalvdf",dir,ionspec)
} else { if (param eq "eedf") {
   filename_base=sprintf("%selectrontotalvdf",dir)
} else { if (param eq "iedf") {
   filename_base=sprintf("%s%stotalvdf",dir,ionspec)
} else {
   print "Unrecognized parameter ".param
   exit
}}}}}}}}}}}}
print "from ".filename_base

do for [i=0:int(maxnum)] {
   filename = sprintf("%s%d.dat",filename_base,i)
   t = system(sprintf("head -%d %s | grep -oE 'SOLUTIONTIME=.{0,20}' | cut -d'=' -f2",numhead,filename)) + 0
   set title leglab
   plot filename\
      using xcol:ycol\
      every ::numhead\
      title sprintf("t = %.3e s",t)
}
