
reset

output_format = "screen"

# ================================
# Terminal and Output Settings
# ================================
if (output_format eq "png") {
    set terminal pngcairo size 1400,600 enhanced font "Arial,16" linewidth 2
    set output "SS_energy_vs_distance.png"
} else { if (output_format eq "pdf") {
    set terminal pdfcairo enhanced color font "Arial,14" linewidth 1.5
    set output "SS_energy_vs_distance.pdf"
} else {
    set terminal wxt 1 size 1400,600 enhanced font "Arial,14"
    set output
}}

input_data_file = "phi1_chi_eBest_SS.dat"

set key off
set border linewidth 2

set yrange [-20:20]     # φ₁ range
set xrange [-230:130]   # χ range (φ₁ + φ₂ | φ₁ - φ₂)

set xtics -200,40,120 font "Arial,11"
set ytics -20,10,20 font "Arial,11" offset 2.5,0
set mxtics 2  # Minor ticks
set mytics 2

set grid xtics ytics mxtics mytics linewidth 1 linecolor rgb "gray70"

set size ratio 1.0/6

set view map

zeta_maximum_value = 1.39687500
energy_reference = -6227.1749

ylabel_definition = ' "φ₁ (°)" font "Arial,12" offset 1.5,0.5'
xlabel_properties = ' "χ (°)" font "Arial,12" offset 0,0.5'
title_properties = 'font "Arial,12" offset 0,-0.6'

colorbox_size = 'size 0.02,0.25'
colorbox_origin_1 = 'origin 0.85,0.70'
colorbox_origin_2 = 'origin 0.85,0.40'
colorbox_origin_3 = 'origin 0.85,0.10'

shift_zeta_for_continuity(zeta) = (zeta < zeta_maximum_value/2) ? zeta : zeta - zeta_maximum_value
shift_phi_to_symmetric_range(phi) = (phi < 180) ? phi : phi - 360

set lmargin at screen 0.05
set rmargin at screen 0.95

set multiplot 

# ===================================================================
# TOP PANEL: RELATIVE BINDING ENERGY SURFACE
# ===================================================================
set tmargin at screen 0.95
set bmargin at screen 0.70

set format x ""
set xlabel ""

set ylabel @ylabel_definition

set title "Relative Binding Energy (eV)" @title_properties

set cbrange [*:*]
set colorbox user @colorbox_origin_1 @colorbox_size

splot input_data_file using 1:2:($5-2*energy_reference) with pm3d notitle

# ===================================================================
# MIDDLE PANEL: STRUCTURAL PARAMETER DISTRIBUTION
# ===================================================================
set tmargin at screen 0.65
set bmargin at screen 0.40

set format x ""
set xlabel ""
set ylabel @ylabel_definition

set title "Optimum Distance (Å)" @title_properties

set palette defined (\
    8.1  '#440154', \
    8.4  '#31688e', \
    8.6  '#21918c', \
    8.8  '#35b779', \
    9.1  '#73d055', \
    9.4  '#bade28', \
    9.7  '#fde725' \
)

set colorbox user @colorbox_origin_2 @colorbox_size

splot input_data_file using 1:2:4 with pm3d notitle

# ===================================================================
# BOTTOM PANEL: ZETA DIHEDRAL ANGLE DISTRIBUTION
# ===================================================================
set tmargin at screen 0.35
set bmargin at screen 0.10

set format x "%.0f"
set xlabel @xlabel_properties
set ylabel @ylabel_definition

set title "ζ Height Difference (Å)" @title_properties

set palette defined (\
    0 '#FF0000', \
    0.5 '#008800', \
    0.9 '#0000FF', \
    1.0 '#FFFFFF', \
    1.1 '#888800', \
    1.5 '#FF00FF', \
    2 '#FF0000'\
)

set colorbox user @colorbox_origin_3 @colorbox_size

# Zeta angle surface with continuity correction
splot input_data_file using 1:2:(shift_zeta_for_continuity($3)) with pm3d notitle

unset multiplot
# ===================================================================