reset

set terminal wxt size 600,800

set grid

set yrange [0:15]
set ytics 0,3,15

cumulative_sum(x)=(a=a+x,a)
set key off

set multiplot

set lmargin at screen 0.1
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.5


set format x ""

a=0; plot "RDF_all.dat" u 1:(cumulative_sum($2+$3+$4+$5)) w l lt rgb '#333333' lw 1.5
set format y ""
unset grid
a=0; plot "RDF_all.dat" u 1:(cumulative_sum(   $2)) w l lt rgb '#FF0000'
a=0; plot "RDF_all.dat" u 1:(cumulative_sum(   $3)) w l lt rgb '#008800'
a=0; plot "RDF_all.dat" u 1:(cumulative_sum(   $4)) w l lt rgb '#0000FF'
a=0; plot "RDF_all.dat" u 1:(cumulative_sum(   $5)) w l lt rgb '#000000'

set lmargin at screen 0.1
set rmargin at screen 0.9
set tmargin at screen 0.5
set bmargin at screen 0.1


set yrange [0:12]
set ytics 0,3,9

set format x "%g"
set format y "%g"

set grid

plot\
	"RDF_all.dat" u 1:(10*($2+$3+$4+$5)) w l lt rgb '#333333' lw 1.5 ,\
	"RDF_all.dat" u 1:(10*(   $2)) w l lt rgb '#FF0000',\
	"RDF_all.dat" u 1:(10*(   $3)) w l lt rgb '#008800',\
	"RDF_all.dat" u 1:(10*(   $4)) w l lt rgb '#0000FF',\
	"RDF_all.dat" u 1:(10*(   $5)) w l lt rgb '#000000',\
	
unset multiplot
