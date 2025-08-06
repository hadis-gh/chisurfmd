set key off

set terminal wxt 1 size 600,800

set xrange [0:1.5*18]


do for[i=1:17]{
	set arrow i from (1.5*i) , graph 0 to (1.5*i) , graph 1 nohead 
}

set palette define (0 '#FF0000',1 '#0000FF')

set multiplot layout 1,2


plot \
	for[i=0:17] sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100   )%360 , 0) u ($5-2*(-6227.1749)):($3+i*1.5):4 w lp pt  9 ps 0.5  palette


plot \
	for [j=0:19:20] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_DD.dat", (i*100+0)%360 , j) u ($5-2*(-6227.1749)):($3+i*1.5):4 w lp pt 7 ps 0.5  palette



#phi1=10
#plot \
#	for [j=0:19:5] for[i=0:17]  sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100+phi1)%360 , j) u ($3+i*1.5):($5-2*(-6227.1749)):4 w lp pt 7 ps 0.5  palette,\


unset multiplot
