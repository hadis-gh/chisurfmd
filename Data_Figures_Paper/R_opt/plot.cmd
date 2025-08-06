
set  grid

set xrange [7.5:10.5]
set yrange [-1.5:0]


set xlabel "distance (Ang)"


set key top left

plot \
	"E_Ropt_minmax_SS.dat" u 1:($2-2*(-6227.1749)) w l lt rgb '#000000' title "SS", "" u 1:($3-2*(-6227.1749)) w l lt rgb '#000000' notitle ,\
	"E_Ropt_minmax_SD.dat" u 1:($2-2*(-6227.1749)) w l lt rgb '#FF0000' title "SD", "" u 1:($3-2*(-6227.1749)) w l lt rgb '#FF0000' notitle ,\
	"E_Ropt_minmax_DS.dat" u 1:($2-2*(-6227.1749)) w l lt rgb '#008800' title "DS", "" u 1:($3-2*(-6227.1749)) w l lt rgb '#008800' notitle ,\
	"E_Ropt_minmax_DD.dat" u 1:($2-2*(-6227.1749)) w l lt rgb '#0000FF' title "DD", "" u 1:($3-2*(-6227.1749)) w l lt rgb '#0000FF' notitle ,\
