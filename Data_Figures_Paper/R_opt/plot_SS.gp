reset

# =============================================================================
# CONFIGURATION PARAMETERS - Modify these for different configurations
# =============================================================================

# Reference energy for normalization
ref_energy = -6227.1749

# Interesting configurations to highlight (phi1, phi2, zeta)
config1_phi1 = 0
config1_phi2 = 0  
config1_zeta = 0.000000
config1_color = "#2E3440"
config1_label = sprintf("φ₁: %d, χ: %d, ζ: %.1f", config1_phi1, config1_phi1 + config1_phi2, config1_zeta)

config2_phi1 = 0
config2_phi2 = 0
config2_zeta = 1.210625
config2_color = "#BF616A"
config2_label = sprintf("φ₁: %d, χ: %d, ζ: %.1f", config2_phi1, config2_phi1 + config2_phi2, config2_zeta)

config3_phi1 = 48
config3_phi2 = 5
config3_zeta = 1.117500
config3_color = "#5E81AC"
config3_label = sprintf("φ₁: %d, χ: %d, ζ: %.1f", config3_phi1, config3_phi1 + config3_phi2, config3_zeta)

config4_phi1 = 196
config4_phi2 = 16
config4_zeta = 0.000000
config4_color = "#A3BE8C"
config4_label = sprintf("φ₁: %d, χ: %d, ζ: %.1f", config4_phi1, config4_phi1 + config4_phi2, config4_zeta)

# Plot appearance settings for publication
gray_points_color = "#CCCCCC"
point_size_gray = 0.8
point_size_config = 1.2
line_width = 2
sampling_interval = 10

# =============================================================================
# PLOT SETUP - Publication Quality Settings
# =============================================================================

# Set output format for publication (uncomment for high-quality output)
# set terminal pngcairo enhanced size 800,600 font "Arial,14"
# set output "SS_interaction_plot.png"

# For screen display with high quality
set terminal wxt enhanced size 900,600 font "Arial,12"

# Enable legend/key display for publication
set key outside right top font "Arial,11" spacing 1.2

# Define axis ranges
set yrange [-1.5:0]
set xrange [8:11]

# Set axis labels with larger fonts
set xlabel "Distance (Å)" font "Arial,14" enhanced
set ylabel "Energy (eV)" font "Arial,14" enhanced

# Set title with publication font
set title "Single-Single (SS) Interaction" font "Arial,16" enhanced

# Improve tick marks and grid
set tics font "Arial,12"
set mxtics 5
set mytics 5
set grid xtics ytics mxtics mytics linewidth 0.5 linecolor rgb "#E0E0E0"

# Set border width
set border linewidth 1.5

# Improve legend appearance
set key box linewidth 1.5 opaque

# =============================================================================
# PLOT: Single-Single (SS) data - Publication Quality
# =============================================================================

plot \
    "../E_Ropt_SS.dat" using 4:($5-2*ref_energy) every sampling_interval \
        pointtype 7 pointsize point_size_gray linecolor rgb gray_points_color \
        title "Optimum distances (all configurations)", \
    "../E_all_SS.dat" using (($1==config1_phi1 && $2==config1_phi2 && $3==config1_zeta) ? $4 : 1/0):(($1==config1_phi1 && $2==config1_phi2 && $3==config1_zeta) ? $5-2*ref_energy : 1/0) \
        with linespoints pointtype 6 pointsize point_size_config linewidth line_width \
        linecolor rgb config1_color title config1_label, \
    "../E_all_SS.dat" using (($1==config2_phi1 && $2==config2_phi2 && $3==config2_zeta) ? $4 : 1/0):(($1==config2_phi1 && $2==config2_phi2 && $3==config2_zeta) ? $5-2*ref_energy : 1/0) \
        with linespoints pointtype 8 pointsize point_size_config linewidth line_width \
        linecolor rgb config2_color title config2_label, \
    "../E_all_SS.dat" using (($1==config3_phi1 && $2==config3_phi2 && $3==config3_zeta) ? $4 : 1/0):(($1==config3_phi1 && $2==config3_phi2 && $3==config3_zeta) ? $5-2*ref_energy : 1/0) \
        with linespoints pointtype 10 pointsize point_size_config linewidth line_width \
        linecolor rgb config3_color title config3_label, \
    "../E_all_SS.dat" using (($1==config4_phi1 && $2==config4_phi2 && $3==config4_zeta) ? $4 : 1/0):(($1==config4_phi1 && $2==config4_phi2 && $3==config4_zeta) ? $5-2*ref_energy : 1/0) \
        with linespoints pointtype 12 pointsize point_size_config linewidth line_width \
        linecolor rgb config4_color title config4_label

# Uncomment these lines to save as high-quality PNG for publication:
# set output
# set terminal wxt
