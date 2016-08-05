set terminal x11
set output
plot "ARGONOutput_pField10.dat"\
   using 1:4\
   every ::3\
   with lines\
   title "argon v_x"
set xlabel "Position (m)"
set ylabel "Velocity (m/s)"
replot
pause -1 "Hit return to quit"
