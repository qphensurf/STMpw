This directory includes some useful utilities and scripts which can be used to process and plot the results:
* Imagen_Bardeen_gnu.f90 and Cond_gnu.f90: fortran programs to process raw data to be plotted using gnuplot.
They can be compiled by doing **make all** in the main directory to generate Imagen_Bardeen_gnu.out and Cond_gnu.out. They are a prerequisite for the provided scripts.
* topo_V_cc.sh (topo_V_ch.sh): script to plot constant current (height) topographic images in the Tersoff-Hamman approximation.
* topo_V_Bardeen_cc.sh (topo_V_Bardeen_ch.sh): same in the Bardeeen approximation.
* cond_V_cc.sh (cond_V_ch.sh): script to plot conductance maps in the constant current (height) mode in the Tersoff-Hamman approximation.
* Imagen_Bardeen_gwy.f90: fortran program to process raw data into a format readable by Gwyddion (https://gwyddion.net)
* dat2siesta.f90: fortran program to convert .dat files into .siesta files readable by WSxM (https://www.wsxm.eu/)
* lut2pal.c: C program to convert a .lut file from the Lut directory of WsXM into .pal to be used as a palette in gnuplot.
* input.STMpw_*: several input.STMpw sample files.
* stm_gui.py: graphical inteface to create STM images and line profiles.
