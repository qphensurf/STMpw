#!/bin/bash
# Script to plot conductance maps in constant current mode

# Variables:

#  Bardeen or TH
   bardeen="F"
   name="Topography_cc"
#   bardeen="F"
#   name="Topography_TH"

#  Voltages to use
   voltages="V_0.500"

#  Currents to use
   currents="1e-9"

#  Repetitions in x and y
   nx=2
   ny=2

#  Values from POSCAR
   a1_x=1.0000000000000000
   a1_y=0.0000000000000000
   a2_x=0.0000000000000000
   a2_y=1.0000000000000000

# End of variables

for j in $voltages
do
cd $j
for i in $currents
do
cat > topo.gnu << !
$bardeen
$nx
$ny
2
$i
0
!
#Imagen_gnu.out < topo.gnu
#Imagen_gnu_Bardeen.out < topo.gnu
Imagen_Bardeen_gnu.out < topo.gnu
cat > stm.gnu << !
set pm3d map
set size ratio -1
set term pngcairo
set output "${name}_${j}_${i}.png"

ex_x=$a1_x; ex_y=$a1_y
ey_x=$a2_x; ey_y=$a2_y

ex2=sqrt(ex_x**2+ex_y**2)
ex_x=ex_x/ex2; ex_y=ex_y/ex2

ey2=sqrt(ey_x**2+ey_y**2)
ey_x=ey_x/ey2; ey_y=ey_y/ey2

set size ratio -1
splot "Topography.gnu" u (ex_x*\$1+ey_x*\$2):(ex_y*\$1+ey_y*\$2):(\$3) notitle
!
gnuplot stm.gnu
done
rm -f stm.gnu Topography.gnu topo.gnu
cd ..
done
