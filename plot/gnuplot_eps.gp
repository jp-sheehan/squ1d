set terminal x11
set output
plot "electronParticles10.dat"\
   using 1:3\
   every ::3\
   with dots\
   title "electrons v_x"
set pointsize 5
set xlabel "Position (m)"
set ylabel "Velocity (m/s)"
replot
pause -1 "Hit return to quit"
