set palette defined ( 0 0 0.5333 0, 0.3333 0 0 1, 0.3333 0 0 0,\
     0.5 0.5333 0 0.5333, 0.6667 1 0 0, 1 1 0.9765 0 )
R = 3.3
set cbrange [-1 : 2]
set xrange [-20 : 140]
set yrange [-20 : 140]
set size square
set key off
plot\
	"bcluster_2.dat" \
	    u 1:2:(R):4 palette with circles fs solid  noborder ,\
	""  u 1:2:(R):($5-1.01)/2 palette lw 1 with circles  ,\
	""  u 1:2:(cos($3/180*3.141592653)*R):(sin($3/180*3.141592653)*R):($5-1.01)/2 with vectors nohead lw 1 palette,\

\
	""  u 1:($2-120):(R):4 palette with circles fs solid  noborder ,\
	""  u 1:($2-120):(R):($5-1.01)/2 palette lw 9 with circles  ,\
	""  u 1:($2-120):(cos($3/180*3.141592653)*R):(sin($3/180*3.141592653)*R):($5-1.01)/2 with vectors nohead lw 9 palette,\
\
	""  u 1:($2+120):(R):4 palette with circles fs solid  noborder ,\
	""  u 1:($2+120):(R):($5-1.01)/2 palette lw 9 with circles  ,\
	""  u 1:($2+120):(cos($3/180*3.141592653)*R):(sin($3/180*3.141592653)*R):($5-1.01)/2 with vectors nohead lw 9 palette,\
\
	""  u ($1-120):($2):(R):4 palette with circles fs solid  noborder ,\
	""  u ($1-120):($2):(R):($5-1.01)/2 palette lw 9 with circles  ,\
	""  u ($1-120):($2):(cos($3/180*3.141592653)*R):(sin($3/180*3.141592653)*R):($5-1.01)/2 with vectors nohead lw 9 palette,\
\
	""  u ($1+120):($2):(R):4 palette with circles fs solid  noborder ,\
	""  u ($1+120):($2):(R):($5-1.01)/2 palette lw 9 with circles  ,\
	""  u ($1+120):($2):(cos($3/180*3.141592653)*R):(sin($3/180*3.141592653)*R):($5-1.01)/2 with vectors nohead lw 9 palette,\
\

