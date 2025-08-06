
reset

 set key off


 set terminal wxt size  800,400

set yrange [-1.5:0]
set xrange [8:11]


set multiplot layout 1,2

plot \
	"../E_Ropt_SS.dat" u 4:($5-2*(-6227.1749)) every 10 pt 7 ps 0.5 lc rgb "#FA000000" ,\
	"<awk '{if($1 ==  0) if($2== 0) if($3==0.000000) print }' ../E_Ropt_SS.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#008800",\
	"<awk '{if($1 ==  0) if($2== 0) if($3==1.210625) print }' ../E_Ropt_SS.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#FF0000",\
	"<awk '{if($1 == 48) if($2== 5) if($3==1.117500) print }' ../E_Ropt_SS.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#0000FF",\
	"<awk '{if($1 ==196) if($2==16) if($3==0.000000) print }' ../E_Ropt_SS.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#FF0000",\

plot \
	"../E_Ropt_DD.dat" u 4:($5-2*(-6227.1749)) every 10 pt 7 ps 0.5 lc rgb "#FA000000" ,\
	"<awk '{if($1 ==  0) if($2== 0) if($3==0.000000) print }' ../E_Ropt_DD.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#008800",\
	"<awk '{if($1 ==  0) if($2== 0) if($3==1.210625) print }' ../E_Ropt_DD.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#FF0000",\
	"<awk '{if($1 == 48) if($2== 5) if($3==1.117500) print }' ../E_Ropt_DD.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#0000FF",\
	"<awk '{if($1 ==196) if($2==16) if($3==0.000000) print }' ../E_Ropt_DD.dat" u 4:($5-2*(-6227.1749)) w lp pt 5 ps 0.3 lc rgb "#FF0000",\


unset multiplot
