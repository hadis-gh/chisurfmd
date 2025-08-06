set key off

set terminal wxt 1 size 800,800





set palette define (0 '#FF0000',1 '#0000FF')

unset colorbox

set multiplot 

#############################################

set lmargin at screen 0.1
set rmargin at screen 0.5
set tmargin at screen 0.9
set bmargin at screen 0.5

do for[i=1:17]{
	set arrow i from (1.5*i) , graph 0 to (1.5*i) , graph 1 nohead 
}

plot \
	for [j=0:19:5] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u ($3+i*1.5):($5-2*(-6227.1749)):4 w lp pt 7 ps 0.5  palette,\

#############################################
set lmargin at screen 0.5
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.5

do for[i=1:17]{
	unset arrow 
}

plot \
	for [j=0:19:20] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u 4:($5-2*(-6227.1749)):($3+i*1.5) w lp pt 7 ps 0.5  palette,\

#############################################
set lmargin at screen 0.5
set rmargin at screen 0.9
set tmargin at screen 0.5
set bmargin at screen 0.1

do for[i=1:17]{
	unset arrow 
}

plot \
	for [j=0:19:20] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u 4:($3+i*1.5):($5-2*(-6227.1749)) w lp pt 7 ps 0.5  palette,\




unset multiplot
