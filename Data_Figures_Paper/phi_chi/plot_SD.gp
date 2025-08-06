reset

set terminal wxt 1 size 1200,400

FILE="phi1_chi_eBest_SD.dat"

set key off

set xrange [-230:130]
set yrange [-20:20]


set format x ""

set xtics -230,20, 130
set ytics -30,20,30

set grid

set size ratio 1.0/6

set view map

zmax=1.39687500
shift_z(x)= (x < zmax/2) ? x : x-zmax
shift_phi(x)= (x < 180) ? x : x-360


set lmargin at screen 0.01
set rmargin at screen 0.99

set multiplot 

set tmargin at screen 0.97
set bmargin at screen 0.67

splot \
	FILE u 1:2:($5-2*(-6227.1749)) with pm3d,\
	
	
set tmargin at screen 0.66
set bmargin at screen 0.36

splot \
	FILE u 1:2:4 with pm3d


set tmargin at screen 0.35
set bmargin at screen 0.05

set palette define (0 '#FF0000' , 0.5 '#008800', 0.9 '#0000FF',1.0 '#FFFFFF',1.1 '#888800',1.5 '#FF00FF', 2 '#FF0000')
splot \
	FILE u 1:2:(shift_z($3)) with pm3d
	
	
	unset multiplot

