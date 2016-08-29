plot 'elec1Particles'.i.'.dat' using 3:(sqrt($4*$4+$5*$5)) title sprintf("speed: i=%i",i)
if (i<n) reread
