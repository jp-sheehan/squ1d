set terminal x11
set output
plot "electronOutput_pField10.dat"\
   using 1:15\
   every ::3\
   with lines\
   title "phi"
set xlabel "Position (m)"
set ylabel "Electric Potential (V)"
replot
pause -1 "Hit return to quit"
