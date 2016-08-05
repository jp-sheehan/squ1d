set terminal x11
set output
plot "electronOutput_pField10.dat"\
   using 1:3\
   every ::3\
   with lines\
   title "electrons"
set xlabel "Position (m)"
set ylabel "Number density (m^-3)"
replot
pause -1 "Hit return to quit"
