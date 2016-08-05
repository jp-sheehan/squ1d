set terminal x11
set output
plot "electronOutput_pField10.dat"\
   using 1:16\
   every ::3\
   with lines\
   title "B_z"
set xlabel "Position (m)"
set ylabel "Magnetic Field (T)"
replot
pause -1 "Hit return to quit"
