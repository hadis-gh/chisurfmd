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
config1_color = "#008800"
config1_label = "Reference (0,0,0.0)"

config2_phi1 = 0
config2_phi2 = 0
config2_zeta = 1.210625
config2_color = "#FF0000"
config2_label = "Config (0,0,1.21)"

config3_phi1 = 48
config3_phi2 = 5
config3_zeta = 1.117500
config3_color = "#0000FF"
config3_label = "Config (48,5,1.12)"

config4_phi1 = 196
config4_phi2 = 16
config4_zeta = 0.000000
config4_color = "#FF0000"
config4_label = "Config (196,16,0.0)"

# Plot appearance settings
gray_points_color = "#FA000000"
point_size_gray = 0.5
point_size_config = 0.3
sampling_interval = 10

# =============================================================================
# PLOT SETUP
# =============================================================================

# Enable legend/key display for publication
set key outside right top

# Set terminal to X11 window with specified dimensions
set terminal wxt size 800,600

# Define axis ranges
set yrange [-1.5:0]
set xrange [8:11]

# Set axis labels
set xlabel "Distance (Å)"
set ylabel "Energy (eV)"

# Set title
set title "Double-Double (DD) Configurations"

# =============================================================================
# PLOT: Double-Double (DD) data
# =============================================================================

plot \
    "../E_Ropt_DD.dat" using 4:($5-2*ref_energy) every sampling_interval pointtype 7 pointsize point_size_gray linecolor rgb gray_points_color title "Optimum distances (all configs)", \
    "../E_all_DD.dat" using (($1==config1_phi1 && $2==config1_phi2 && $3==config1_zeta) ? $4 : 1/0):(($1==config1_phi1 && $2==config1_phi2 && $3==config1_zeta) ? $5-2*ref_energy : 1/0) with linespoints pointtype 5 pointsize point_size_config linecolor rgb config1_color title config1_label, \
    "../E_all_DD.dat" using (($1==config2_phi1 && $2==config2_phi2 && $3==config2_zeta) ? $4 : 1/0):(($1==config2_phi1 && $2==config2_phi2 && $3==config2_zeta) ? $5-2*ref_energy : 1/0) with linespoints pointtype 5 pointsize point_size_config linecolor rgb config2_color title config2_label, \
    "../E_all_DD.dat" using (($1==config3_phi1 && $2==config3_phi2 && $3==config3_zeta) ? $4 : 1/0):(($1==config3_phi1 && $2==config3_phi2 && $3==config3_zeta) ? $5-2*ref_energy : 1/0) with linespoints pointtype 5 pointsize point_size_config linecolor rgb config3_color title config3_label, \
    "../E_all_DD.dat" using (($1==config4_phi1 && $2==config4_phi2 && $3==config4_zeta) ? $4 : 1/0):(($1==config4_phi1 && $2==config4_phi2 && $3==config4_zeta) ? $5-2*ref_energy : 1/0) with linespoints pointtype 5 pointsize point_size_config linecolor rgb config4_color title config4_label
