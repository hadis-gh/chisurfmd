reset

output_format = "screen"

# ================================
# Terminal and Output Settings
# ================================
if (output_format eq "png") {
    set terminal pngcairo size 1200,1200 enhanced font "Arial,20"
    set output "DS_energy_vs_distance.png"
} else { if (output_format eq "eps") {
    set terminal postscript eps size 8,8 enhanced color font "Arial,20"
    set output "DS_energy_vs_distance.eps"
} else {
    set terminal wxt 1 size 900,900 enhanced font "Arial,18"
    set output
}}

set title "DS Interaction: Energy Distance" font "Arial,22" offset 0,-0.5

# ================================
# Data File Configuration
# ================================
ref_energy = -6227.1749

input_data_file = "E_Ropt_minmax_DS.dat"
gray_data_file = "../E_Ropt_DS.dat"
all_data_file = "../E_all_DS.dat"

# ================================
# Plot Layout and Styling
# ================================
set size ratio 1.0
set key inside right bottom box opaque font "Arial,22" spacing 1.2 maxcols 1 width -10 samplen .5 offset 0,.8
set border linewidth 2.0
set tics scale 1.2

set grid xtics ytics mxtics mytics linewidth 1.0 linecolor rgb "gray70" dashtype 2

# ================================
# Axis Configuration
# ================================
x_min = 7.5
x_max = 11.0
set xrange [x_min:x_max]

y_min = -1.5
y_max = 0.0
set yrange [y_min:y_max]

set xlabel "Distance (Å)" offset 0,0.5 font "Arial,22"
set ylabel "Binding Energy (eV)" offset 0,0 font "Arial,22"

set xtics x_min,0.5,x_max font "Arial,20"
set ytics auto font "Arial,20"
set mxtics 2
set mytics 2

# ================================
# Highlighted Configurations
# ================================
array config_phi1[4] = [233,    50.0,        0.0,        131]
array config_chi[4] =  [229,    40.0,        0.0,       0]
array config_zeta[4] = [0.093125,   1.024375,   0.0,    0.558750]

array config_phi2[4]
do for [i=1:4] {
    config_phi2[i] = config_phi1[i] - config_chi[i]
}

array config_colors[4] = ["#3b4cc0", "#4daf4a", "#ff7f00", "#b40426"]
array config_labels[4]

do for [i=1:4] {
    config_labels[i] = sprintf("φ₁=%.0f°, χ=%.0f°, ζ=%.2f°", \
                              config_phi1[i], config_chi[i], config_zeta[i])
}

# ================================
# All data file for background
# ================================
bg_point_size = 0.4
bg_sampling = 10
bg_color = "#FA000000"

highlight_point_size = 0.8
highlight_line_width = 2.0

# ================================
# Main Plotting Commands
# ================================
plot \
    gray_data_file u 4:($5 - 2*ref_energy) every bg_sampling \
        with points pointtype 7 pointsize bg_point_size \
        linecolor rgb bg_color \
        title "All configurations", \
    for [i=1:4] sprintf('<awk ''$1==%f && $2==%f && $3==%f {print}'' %s', \
                       config_phi1[i], config_phi2[i], config_zeta[i], all_data_file) \
        u 4:($5 - 2*ref_energy) \
        with linespoints pointtype 7 pointsize highlight_point_size \
        linewidth highlight_line_width linecolor rgb config_colors[i] \
        title config_labels[i]
