reset
set terminal wxt 1 size 1400,600 enhanced font "Arial,12"

input_data_file = "phi1_chi_eBest_SS.dat"

set key off
set border linewidth 2

set yrange [-20:20]
set xrange [-230:130]

set xtics -200,40,120 font "Arial,11"
set ytics -20,10,20 font "Arial,11" offset 2.5,0
set mxtics 2
set mytics 2

set grid xtics ytics mxtics mytics linewidth 1 linecolor rgb "gray70"

set size ratio 1.0/6

set view map

zeta_maximum_value = 1.39687500
energy_reference = -6227.1749

ylabel_definition = ' "φ₁ (°)" font "Arial,14" offset 1.5,0.5'
xlabel_properties = ' "χ (°)" font "Arial,14" offset 0,0.5'
title_properties = 'font "Arial,14" offset 0,-0.6'

colorbox_size = 'size 0.02,0.28'
colorbox_origin_1 = 'origin 0.885,0.69'
colorbox_origin_2 = 'origin 0.885,0.39'
colorbox_origin_3 = 'origin 0.885,0.09'

vertical_title_font = 'font "Arial,12"'
energy_title_pos =     'at screen 0.944, 0.84'
distance_title_pos =   'at screen 0.944, 0.54' 
zeta_title_pos =       'at screen 0.944, 0.24'

shift_zeta_for_continuity(zeta) = (zeta < zeta_maximum_value/2) ? zeta : zeta - zeta_maximum_value
shift_phi_to_symmetric_range(phi) = (phi < 180) ? phi : phi - 360

set lmargin at screen 0.08
set rmargin at screen 0.95

set multiplot 

# ===================================================================
# TOP PANEL: RELATIVE BINDING ENERGY SURFACE
# ===================================================================
set tmargin at screen 0.97
set bmargin at screen 0.69

set format x ""
set xlabel ""

set ylabel @ylabel_definition

unset title

set cbrange [*:*]
set colorbox user @colorbox_origin_1 @colorbox_size

unset label
set label 1 "Binding Energy (J)" @energy_title_pos @vertical_title_font tc rgb "black" center rotate by 90

splot input_data_file using 1:2:($5-2*energy_reference) with pm3d notitle

# ===================================================================
# MIDDLE PANEL: STRUCTURAL PARAMETER DISTRIBUTION
# ===================================================================
set tmargin at screen 0.67
set bmargin at screen 0.39

set format x ""
set xlabel ""
set ylabel @ylabel_definition

unset title

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

unset label
set label 2 "Optimum Distance (Å)" @distance_title_pos @vertical_title_font tc rgb "black" center rotate by 90

splot input_data_file using 1:2:4 with pm3d notitle

# ===================================================================
# BOTTOM PANEL: ZETA DIHEDRAL ANGLE DISTRIBUTION
# ===================================================================
set tmargin at screen 0.37
set bmargin at screen 0.09

set format x "%.0f"
set xlabel @xlabel_properties
set ylabel @ylabel_definition

unset title

set palette define (0 '#FF0000', 0.5 '#008800', 0.9 '#0000FF', 1.0 '#FFFFFF', 1.1 '#888800', 1.5 '#FF00FF', 2 '#FF0000')

set colorbox user @colorbox_origin_3 @colorbox_size

unset label
set label 3 "ζ Height Difference (Å)" @zeta_title_pos @vertical_title_font tc rgb "black" center rotate by 90

splot input_data_file using 1:2:(shift_zeta_for_continuity($3)) with pm3d notitle

unset multiplot