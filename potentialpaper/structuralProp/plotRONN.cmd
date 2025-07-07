### horizontal histogram
reset session

set key off

set style fill transparent solid 0.5 noborder
set multiplot 

set grid

####################################### CHECK NUMBER OF POINTS

stats "<grep -e 'SS' data_all.dat"   nooutput; if (!exists("STATS_records")) { NumSS=0 ;} else {NumSS=STATS_records}
stats "<grep -e 'SD' data_all.dat"   nooutput; if (!exists("STATS_records")) { NumSD=0 ;} else {NumSD=STATS_records}
stats "<grep -e 'DS' data_all.dat"   nooutput; if (!exists("STATS_records")) { NumDS=0 ;} else {NumDS=STATS_records}
stats "<grep -e 'DD' data_all.dat"   nooutput; if (!exists("STATS_records")) { NumDD=0 ;} else {NumDD=STATS_records}

print "SS-interaction= ", NumSS
print "SD-interaction= ", NumSD
print "DS-interaction= ", NumDS
print "DD-interaction= ", NumDD

####################################### PLOT SCATTER

set lmargin at screen 0.1
set rmargin at screen 0.65
set tmargin at screen 0.65
set bmargin at screen 0.1

set xrange [-230:130]
set yrange [0:1.5]

plot \
	"<grep -e 'SS' data_all.dat" u 10:6 pt 7 ps 0.1 lc '#FF0000',\
	"<grep -e 'SD' data_all.dat" u 10:6 pt 7 ps 0.1 lc '#0000FF',\
	"<grep -e 'DS' data_all.dat" u 10:6 pt 7 ps 0.1 lc '#008800',\
	"<grep -e 'DD' data_all.dat" u 10:6 pt 7 ps 0.1 lc '#000000',\

####################################### PLOT HISTOGRAM A



set x2range [-230:130]
set yrange [0:*]

set ytics textcolor rgb "#FFFFFF00"

unset xtics
set x2tics
set grid x2tics

set lmargin at screen 0.1
set rmargin at screen 0.65
set tmargin at screen 0.9
set bmargin at screen 0.65

binwidth = 1
bin(x)   = binwidth * floor(x/binwidth)

if(NumSS!=0){ set table $HistoA_1 ; plot "<grep -e 'SS' data_all.dat" u (bin($10)):(1.0/NumSS) smooth freq ; unset table } else { set table $HistoA_1 ; plot '+' u 1:(0.0) ; unset table }
if(NumSD!=0){ set table $HistoA_2 ; plot "<grep -e 'SD' data_all.dat" u (bin($10)):(1.0/NumSD) smooth freq ; unset table } else { set table $HistoA_2 ; plot '+' u 1:(0.0) ; unset table }
if(NumDS!=0){ set table $HistoA_3 ; plot "<grep -e 'DS' data_all.dat" u (bin($10)):(1.0/NumDS) smooth freq ; unset table } else { set table $HistoA_3 ; plot '+' u 1:(0.0) ; unset table }
if(NumDD!=0){ set table $HistoA_4 ; plot "<grep -e 'DD' data_all.dat" u (bin($10)):(1.0/NumDD) smooth freq ; unset table } else { set table $HistoA_4 ; plot '+' u 1:(0.0) ; unset table }

plot \
	$HistoA_1 u 1:($2/2):(binwidth/2.):($2/2) w boxxy axis x2y1 lc '#FF0000',\
	$HistoA_2 u 1:($2/2):(binwidth/2.):($2/2) w boxxy axis x2y1 lc '#0000FF',\
	$HistoA_3 u 1:($2/2):(binwidth/2.):($2/2) w boxxy axis x2y1 lc '#008800',\
	$HistoA_4 u 1:($2/2):(binwidth/2.):($2/2) w boxxy axis x2y1 lc '#000000'


####################################### PLOT HISTOGRAM Z

set xrange [0:1.1]
set y2range [0:1.5]

set xtics textcolor rgb "#FFFFFF00"

unset ytics
unset x2tics
set y2tics
set grid y2tics

set lmargin at screen 0.65
set rmargin at screen 0.9
set tmargin at screen 0.65
set bmargin at screen 0.1


binwidth = 0.1
bin(x)   = binwidth * floor(x/binwidth)
if(NumSS!=0){ set table $HistoZ_1 ; plot "<grep -e 'SS' data_all.dat" u (bin( $6)):(1.0/NumSS) smooth freq ; unset table } else { set table $HistoZ_1 ; plot '+' u 1:(0.0) ; unset table }
if(NumSD!=0){ set table $HistoZ_2 ; plot "<grep -e 'SD' data_all.dat" u (bin( $6)):(1.0/NumSD) smooth freq ; unset table } else { set table $HistoZ_2 ; plot '+' u 1:(0.0) ; unset table }
if(NumDS!=0){ set table $HistoZ_3 ; plot "<grep -e 'DS' data_all.dat" u (bin( $6)):(1.0/NumDS) smooth freq ; unset table } else { set table $HistoZ_3 ; plot '+' u 1:(0.0) ; unset table }
if(NumDD!=0){ set table $HistoZ_4 ; plot "<grep -e 'DD' data_all.dat" u (bin( $6)):(1.0/NumDD) smooth freq ; unset table } else { set table $HistoZ_4 ; plot '+' u 1:(0.0) ; unset table }

plot \
	$HistoZ_1 u ($2/2):1:($2/2):(binwidth/2.) w boxxy axis x1y2 lc '#FF0000',\
	$HistoZ_2 u ($2/2):1:($2/2):(binwidth/2.) w boxxy axis x1y2 lc '#0000FF',\
	$HistoZ_3 u ($2/2):1:($2/2):(binwidth/2.) w boxxy axis x1y2 lc '#008800',\
	$HistoZ_4 u ($2/2):1:($2/2):(binwidth/2.) w boxxy axis x1y2 lc '#000000'
	


unset multiplot
