reset

set terminal wxt 1 size 1200,400

FILE="phi1_chi_eBest_DD.dat"

set key off



set grid

#set size ratio 1.0/6

set view map

zmax=1.39687500
shift_z(x)= (x < zmax/2) ? x : x-zmax
shift_phi(x)= (x < 180) ? x : x-360


splot \
	FILE u 1:2:($5-2*(-6227.1749)) with pm3d,\


