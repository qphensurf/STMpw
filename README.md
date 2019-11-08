# STMpw

STMpw is a program to perform STM constant current and dI/dV map simulations in the *Tersoff-Hamann* and *Bardeen* approximations. It is a postprocessing utility for DFT planewave codes. At the moment it is interfaced with [VASP](https://www.vasp.at).

The program can perform:

    a) Constant current and constant height STM topographic images
    b) dI/dV maps in constant current and constant height modes
    c) dI/dV curves

## Approximations

* (*TH & Bardeen*) The program imposes that the wavefunctions exponentially decay into vacuum beyond a certain distance **z_s** into vacuum. Starting by a plane parallel to the surface at a distance **z_s** located above the system of interest each of the planewave components of the wavefunctions from the DFT calculation are matched to an exponential function decaying into vacuum. If the surface potential is already zero at **z_s** this procedure is exact.

* (*TH & Bardeen*) Temperature is zero. The bias voltage effect will be much larger.

* (*Bardeen*) The energy-conservation delta's are replaced by Gaussians with a hardwired *sigma* of 0.25 eV.

* (*Bardeen*) The wavefunctions of the tip are also matched to exponential functions from a plane **z_t**.

* (*Bardeen*) In order to make it very fast, the main approximation is to calculate a tip in the same cell as the object under study.

* (*Bardeen*) The program assumes the same workfunction for tip and sample. This approximation is not necessary, but makes life easier and our experience is that, in first approximation, it's more than fine.

## Compilation

The program is compiled using the included Makefile by typing *make* in the main directory. The program has been tested in linux with INTEL and GNU FORTRAN compilers. 

It is required to provide a FFTW implementation. The performance of the program will dramatically depend on the performance of the FFTW library, so it is strongly recommended to use an optimized version for your system. MKL typically offers a good performance.
 
The Makefile file must be adapted to your system by choosing a FORTRAN compiler and a FFTW library.

## Usage

1.  You must run VASP with some requirements:

* Schematically the unit cell must have the form:

             |ttttttttttt|
             |ttttttttttt|
             |    ttt    |
             |     t     |
             |           |
             |           |
             |           |
             |    mmm    |
             |sssssssssss|
             |sssssssssss|
             |sssssssssss|
             
        t: tip atoms (for Bardeen only)
        s: substrate atoms
        m: molecule atoms

	Therefore the molecule must be above the surface and (for Bardeen) the tip must be above the molecule. While the surface unit cell can have any shape, the vertical axis must be perpendicular to the surface, in the following all distances along this axis are called **z**.

* VASP must be run in the **standard** version (the gamma version is not supported at this moment, if you run gamma-only calculations you can just rerun with the standard vasp version including NELM=0 in the INCAR file, this quickly produces a good input for the code).
* The use of symmetry in VASP is not allowed: SYM=0 **must** be included in the INCAR file.
* The calculation must be carefully converged in k-points, especially for dI/dV maps and curves.

	Once the VASP has converged POSCAR, OUTCAR and WAVECAR must be kept and supplied to STMpw. 

2. You must determine the plane **z_s** from which the wavefunction Fourier components will be substituted by exponential functions of **z**. First, you have to look for the last surface layer to be used as reference (**Zsurf**). Then you look at the density vs density behavior to look for the first point where there is an exponential decay. The z value of this point in direct coordinates will be **z_s**. Typically this point is around 2-3 angstroms above the most protruding atom of the molecule. So to speed things up a negative value of **z_s** indicates a distance in angstroms above the most protruding atom (the code will look for the most protruding atom). In addition, if you are using the Bardeen approach the origin of the tip **Ztip** must be supplied. **z_t** is hardwired to 2 angstroms below the most (least?) protruding atom of the tip. At the end you should have **Zsurf < z_s < z_t < Ztip**.

3. For the calculation of dI/dV curves you must determine 
* The range of energies (or bias value) to use (referred to Fermi energy = 0 eV). 
* The number of divisions (or total number of bias values)  on that range
* On which points the curves will be determined. Tipically we select the *x* and *y* coordinates of one atom (in angstroms) and a *z* value several angstroms over the *z* of the atom.

4. You must create an input.STMpw file as follows (where ! means that a comment follows, before ! you have to write a value for the marked variable. Please start from the first line of the file.):

   	   phi   ! workfunction of the surface in eV 
	   n   ! number of voltages to calculate
	   V1...Vn   ! n values of tip-substrate voltages in V (tip to mass)
	   nZ   ! sampling points in z (perpendicular to the surface)
	   Zmax   ! maximum tip-surface distance (from 'Zsurf') in angstroms.
	   Zsurf   ! origin of the surface in direct coordinates
	   z_s   ! as explained above
	   Bardeen    ! T for Bardeen and and F for Tersoff-Hamman
	   Ztip   ! (only if Bardeen=T) origin of the tip in direct coordinates
	   dIdV   ! T or F, whether we calculate the dIdV curve
	   emin emax   ! (if dIdV=T) range of energy for dIdV
	   ndiv   ! (if dIdV=T) number of divisions between emin and emax
	   Np ! (if dIdV=T) number of points to plot the dIdV curve
	   x1, y1, z1 ! (if dIdV=T) coordinates of the first point (in angstroms)
	   ...
	   xNp, yNp, zNp ! (if dIdV=T) coordinates of the nP point (in angstroms)
	   name_POSCAR   ! name of the POSCAR file
	   name_WAVECAR   ! name of the WAVECAR file 
	   mapfile   ! T to read a reciprocal vector and index file or F to generate it
	   name_mapfile   ! (if mapfile=T) name of the reciprocal vector and index file
	   Gamma   ! T or F, whether the k-point sampling contains the Gamma point or not.
	   wsxm   ! T or F, WSxM output?
	   factor   ! (if wsxm=T) multiplying factor for WSxM output files
	   gnuplot ! T or F, plain output to use in gnuplot?
	   cube   ! T or F, cube format output?

	Some examples can be found in Utils.

5. Different output files are produced. For each required voltage a directory with the name **V_voltage** is created. Inside each directory, different files can be found depending on the required output. 
* **WSxM**: output for the [WSxM](http://www.wsxm.es/download.html) program. There is a TH_V_voltage.siesta file for *Tersoff-Hamann* and Bardeen_V_voltage.siesta and TH_tip_V_voltage.siesta files for *Bardeen*. The last file is a STM simulation of the tip in the *Tersoff-Hamann* approximation. All these files can be directly read by WsXM (choose the *.* format for reading it in WSxM) and processed using: Process -> Filter -> Create STM type image...
* **gnuplot**: .dat plain files. They can be plotted, for example, with gnuplot. In the Utils directory there are different programs and scripts to process them.
* **cube**: files in the cube format. At the moment there are just for *Tersoff-Hamann*: TH_V_voltage.cube for STM images and dIdV_TH_V_voltage.cube for dIdV maps. cube is a standard format which can be read by many programs, including the last versions of WSxM.

	**Note**: distances are referred to both surfaces (sample and tip) but the sampling region is only between **z_s** and **z_t**, because it is the 'asymptotic' region.

## Authors
Nicolás Lorente and Roberto Robles based on the Bardeen2 code of Nicolás Lorente.

## License
GNU General Public License v3.0
