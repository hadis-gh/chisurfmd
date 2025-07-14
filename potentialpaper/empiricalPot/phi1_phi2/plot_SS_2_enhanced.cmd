# ============================================================================
# Gnuplot Script for phi1-phi2 Energy Surface Analysis
# -----------------------------------------------------------------------------
# Generates three 2D heatmaps (top row) + and three 1D cross-section plots (bottom row)
# ============================================================================

reset   # Reset all gnuplot settings to default

FILE="phi1_phi2_eBest_SS.dat"   # Input data file

set key off                     # Hide legend/key

# -----------------------------------------------------------------------------
# Axis and grid settings for all panels
# -----------------------------------------------------------------------------
set xrange [-180:180]         # X axis range (phi1)
set xtics -180,60,180         # X axis ticks every 60 degrees
set ytics -180,60,180         # Y axis ticks every 60 degrees
set yrange [-180:180]         # Y axis range (phi2)
set format x ""               # Hide x axis numbers for top panels
set grid                      # Show grid lines
set view map                  # 2D map view
set size square               # Square aspect ratio

set terminal wxt size 1200,600 # Output window size for preview
set multiplot layout 1,3       # Arrange three plots side by side in one row

# -----------------------------------------------------------------------------
# Helper functions for shifting data (for periodicity)
# -----------------------------------------------------------------------------
zmax=1.39687500
shift_z(x)= (x < zmax/2) ? x : x-zmax   # Shift z parameter for periodicity
shift_phi(x)= (x < 180) ? x : x-360     # Shift phi for periodicity

# -----------------------------------------------------------------------------
# Panel 1: Energy Surface Heatmap
# -----------------------------------------------------------------------------
set lmargin at screen 0.1
set rmargin at screen 0.3
set tmargin at screen 0.9
set bmargin at screen 0.4

# Add geometric guide arrows for visual reference
set arrow 1 from -60,-10 to 120,-10 nohead lw 3 lc '#888888' front
set arrow 2 from -40, 10 to 140, 10 nohead lw 3 lc '#888888' front
set arrow 3 from -60,-10 to -40, 10 nohead lw 3 lc '#888888' front
set arrow 4 from 120,-10 to 140, 10 nohead lw 3 lc '#888888' front

# Plot: Energy surface as a heatmap (pm3d)
#   - Uses columns 1 and 2 for phi1 and phi2, column 5 for energy (with offset)
#   - Overlays guide lines for visual reference

splot \
	FILE u 1:2:($5-2*(-6227.1749)) with pm3d,\
	'+' using 1:($1+ 50):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-130):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1+230):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-310):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(-0.5) with lines lc rgb '#CCCCCC' dt 3 lw 1,\

# -----------------------------------------------------------------------------
# Panel 2: Force Field Parameter Heatmap
# -----------------------------------------------------------------------------
set lmargin at screen 0.4
set rmargin at screen 0.6
set tmargin at screen 0.9
set bmargin at screen 0.4

# Plot: Force field parameter as a heatmap (pm3d)
#   - Uses shifted phi1 and phi2, column 4 for parameter

splot FILE u ( shift_phi($1)):( shift_phi($2)):4 with pm3d,\
	'+' using 1:($1+ 50):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-130):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1+230):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-310):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(9) with lines lc rgb '#CCCCCC' dt 3 lw 1,\

# -----------------------------------------------------------------------------
# Panel 3: Z-Parameter Heatmap
# -----------------------------------------------------------------------------
set lmargin at screen 0.7
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.4

set palette define (0 '#FF0000' , 0.5 '#008800', 0.9 '#0000FF',1.0 '#FFFFFF',1.1 '#888800',1.5 '#FF00FF', 2 '#FF0000')          # Custom color palette for Z-parameter

# Plot: Z-parameter as a heatmap (pm3d)
#   - Uses shifted phi1, phi2, and z
splot FILE u ( shift_phi($1)):( shift_phi($2)):(shift_z($3)) with pm3d,\
	'+' using 1:($1+ 50):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-130):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1+230):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	'+' using 1:($1-310):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	for [i=-170:180:20]'+' using 1:(1.0*i):(0) with lines lc rgb '#CCCCCC' dt 3 lw 1,\
	
#########################################################
# End of top row (heatmaps)
#########################################################

unset arrow 1 
unset arrow 2 
unset arrow 3 
unset arrow 4 

# -----------------------------------------------------------------------------
# Bottom row: 1D cross-sections at phi2 = 0
# -----------------------------------------------------------------------------
set size ratio 0.5
set yrange [*:*]
set ytics auto
set format x "%g"

# Panel 4: Energy cross-section
set lmargin at screen 0.1
set rmargin at screen 0.3
set tmargin at screen 0.4
set bmargin at screen 0.2

# Plot: Energy vs phi1 for phi2=0
plot \
	sprintf("<cat %s | awk '{if ($2==0) print }'",FILE) u 1:($5-2*(-6227.1749)) w l lt 1

# Panel 5: Force field parameter cross-section
set lmargin at screen 0.4
set rmargin at screen 0.6
set tmargin at screen 0.4
set bmargin at screen 0.2

# Plot: Parameter vs phi1 for phi2=0
plot \
	sprintf("<cat %s | awk '{if ($2==0) print }'",FILE) u 1:($4) w l lt 1

# Panel 6: Z-parameter cross-section
set lmargin at screen 0.7
set rmargin at screen 0.9
set tmargin at screen 0.4
set bmargin at screen 0.2

# Plot: Z vs phi1 for phi2=0
plot \
	sprintf("<cat %s | awk '{if ($2==0) print }'",FILE) u 1:(shift_z($3)) w l lt 1
unset multiplot
