reset


FILE="phi1_phi2_eBest_DD.dat"

set key off

set xrange [-180:180]
set yrange [-180:180]

set format x ""

set xtics -180,60,180
set ytics -180,60,1800

set grid

set view map

set size square

set terminal wxt size 1200,600
set multiplot layout 1,3

zmax=1.39687500
shift_z(x)= (x < zmax/2) ? x : x-zmax
shift_phi(x)= (x < 180) ? x : x-360


set lmargin at screen 0.1
set rmargin at screen 0.3
set tmargin at screen 0.9
set bmargin at screen 0.4

set arrow 1 from -180,-10 to 180,-10 nohead lw 3 lc '#888888' front
set arrow 2 from -180, 10 to 180, 10 nohead lw 3 lc '#888888' front

splot \
	FILE u 1:2:($5-2*(-6227.1749)) with pm3d,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\


set lmargin at screen 0.4
set rmargin at screen 0.6
set tmargin at screen 0.9
set bmargin at screen 0.4

splot FILE u ( shift_phi($1)):( shift_phi($2)):4 with pm3d,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\

set lmargin at screen 0.7
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.4

set palette define (0 '#FF0000' , 0.5 '#008800', 0.9 '#0000FF',1.0 '#FFFFFF',1.1 '#888800',1.5 '#FF00FF', 2 '#FF0000')
splot FILE u ( shift_phi($1)):( shift_phi($2)):(shift_z($3)) with pm3d,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	
#########################################################
#
#########################################################

unset arrow 1 
unset arrow 2 
unset arrow 3 
unset arrow 4 

set size ratio 0.5
set yrange [*:*]

set ytics auto
set format x "%g"

########################################################
set lmargin at screen 0.1
set rmargin at screen 0.3
set tmargin at screen 0.4
set bmargin at screen 0.2

plot \
	sprintf("<cat %s | awk '{if ($2== 0) print }'",FILE) u 1:($5-2*(-6227.1749)) w l lt 1,\
	sprintf("<cat %s | awk '{if ($2==10) print }'",FILE) u ($1+10):($5-2*(-6227.1749)) w l lt 2 dt 2,\

########################################################
set lmargin at screen 0.4
set rmargin at screen 0.6
set tmargin at screen 0.4
set bmargin at screen 0.2

plot \
	sprintf("<cat %s | awk '{if ($2==0) print }'",FILE) u 1:($4) w l lt 1

########################################################
set lmargin at screen 0.7
set rmargin at screen 0.9
set tmargin at screen 0.4
set bmargin at screen 0.2

plot \
	sprintf("<cat %s | awk '{if ($2==0) print }'",FILE) u 1:(shift_z($3)) w l lt 1,\
	sprintf("<cat %s | awk '{if ($2==5) print }'",FILE) u ($1+1):(shift_z($3)) w l lt 1
unset multiplot
