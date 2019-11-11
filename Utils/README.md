This directory includes some very simple utilities and scripts which can be used to plot the results
using gnuplot. 
* Imagen_Bardeen_gnu.f90 and Cond_gnu.f90: fortran programs to process raw data to be plotted using gnuplot.
They can be compiled by doing **make all** in the main directory to generate Imagen_Bardeen_gnu.out and Cond_gnu.out. They are a prerequisite for the provided scripts.
* topo_V_cc.sh (topo_V_ch.sh): script to plot constant current (height) topographic images in the Tersoff-Hamman approximation.
* topo_V_Bardeen_cc.sh (topo_V_Bardeen_ch.sh): same in the Bardeeen approximation.
* cond_V_cc.sh (cond_V_ch.sh): script to plot conductance maps in the constant current (height) mode in the Tersoff-Hamman approximation.
