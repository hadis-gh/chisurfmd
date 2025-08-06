set key off

set terminal wxt 1 size 800,600

set yrange [0:1.5*18]
set xrange [8:11]

set palette define (0 '#FF0000',1 '#0000FF')

set multiplot layout 1,2


plot \
	for[i=0:17] sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100   )%360 , 0) u 4:($3+i*1.5):($5-2*(-6227.1749)) w lp pt  9 ps 0.5  palette,\
	for[i=0:17] "E_zetaBest_Rbest_SS.dat" u 4:($3+i*1.5):($5-2*(-6227.1749)) w lp pt  9 ps 0.5  palette,\
	for[i=0:17] i*1.5 w l lt -1 

plot \
	for[j=0:19:20] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_DD.dat", (i*100+0)%360 , j) u 4:($3+i*1.5):($5-2*(-6227.1749)) w lp pt 7 ps 0.5  palette,\
	for[i=0:17] "E_zetaBest_Rbest_DD.dat" u 4:($3+i*1.5):($5-2*(-6227.1749)) w lp pt  9 ps 0.5  palette,\
	for[i=0:17] i*1.5 w l lt -1 





unset multiplot
