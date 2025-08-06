
set view 60, 30

set xrange [*:*]
set yrange [*:*]
set zrange [*:*]

set xlabel "X"
set ylabel "Y"
set zlabel "Z"

set zrange [8:10]

j=0

splot \
    for[i=0:17 ]sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u ($3+i*1.5):($5-2*(-6227.1749)):4 with points pt 7 lc rgb "blue" title '3D Points', \
    for[i=0:17 ]sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u ($3+i*1.5):($5-2*(-6227.1749)):(8.0) with points pt 7 lc rgb "gray" notitle, \
    for[i=0:17 ]sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u ($3+i*1.5):(0.0              ):4 with points pt 7 lc rgb "gray" notitle, \
    for[i=0:17 ]sprintf("<awk '{if($1 == %d) if($2== %d)  print }' E_zeta_Rbest_SS.dat", (i*100)%360 , j) u (0.0):($5-2*(-6227.1749)):4 with points pt 7 lc rgb "gray" notitle 
