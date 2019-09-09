STMpw.f90

Simple implementation for Bardeen's transfer hamiltonian
and tunneling extension for the STM

Depending on the paramenters of 
the input file : input.STMpw
this program can perform

      a) Constant current STM images
      b) dI/dV maps 

Approximations:

  This program is to be used as a VASP post-processing 
utility. In order to make it very fast the main approximation is to
calculate a tip in the same cell as the object under study.

  Plus the k-point sampling should be such that no k-points falls
on the Brillouin zone boundary, otherwise the real space integrals
give delta's on the reciprocal vectors that can be satisfied by
combinations of k-points and reciprocal vectors G.

  Third: we have assumed the same workfunction for tip and sample
this approximation is not necessary but makes life easier and our
experience is that in first approximation, it's more than fine.

  Fourth: temperature is zero. The bias voltage effect will be
much larger.

  Fifth: the energy-conservation delta's are replaced by gaussians
following the equation

      \begin{equation}
         \delta (E) \approx 
            \frac{1}{\sqrt(\pi) \sigma} \exp( - (E/\sigma)^2 )
      \end{equation}

\sigma is the variable sigma in this code and is hardwired to 0.25 eV
in declaration.f90

  Sixth: the wavefunctions are exponentially decaying in vacuum, using
the following expression

      \begin{equation}
         \psi(\vect{r}) \sim
            \frac{1}{\sqrt{vol}} \sum_{\vect{G}} A_G (z_s) 
                  \exp(-\sqrt(2\phi+(k+G)^2)|z-z_s|)
                  \exp(i G \cdot r)
      \end{equation}

where G are the reciprocal vectors of the 3-D FFT. We match
this asymptotic form of the wf at a plane located at z_s (in
the code we keep this notation) with workfuntion phi. 
For the tip wavefunction we use the same notation except for
the 2-D FFT are called C_G and the matching distance is
z_t.


Instructions:

First run either VASP (for a and b above) or deltaWF (for c)
where the atomic configuration should be

          |         ttttssss      |
          |        tttttssssm     |
          |       ttttttssssm     |
          |        tttttssssm     |
          |         ttttssss      |
                 
          |<-------- unit cell -->|

      t: 'tip' atoms
      s: 'substrate' atoms
      m: 'molecule' atoms

This code needs the tip higher than the molecule
as is seen in the above scheme if you use the
periodicity of the cell.

Second determine z_s and z_t (see Sixth above) as
the distance to the surface (last layer 's') by
looking at the density versus 'z' behavior: the
first point where it is exponentially decaying near
the molecule (z_s) and the tip (z_t).
If z_s < 0 it is the distance (in A) above the topmost atom
 starting from which we assume an exponential decay.

z_t is hardwired:
z_t is calculated in this code as 2.0 Angs below
the most protruding atom below Ztip.

Hence Zsurf < z_s < z_t < Ztip is compulsory!

ACHTUNG BITTE:
     Distances are then referred to both surfaces (sample and tip)
     but the sampling region is only between z_s and z_t, because it
     is the 'asymptotic' region.

Third fill in the file input.STMpw with the following lines

-----------------------------------------------------------------
                 input.STMpw
-----------------------------------------------------------------


    T ! whether we use Bardeen aprox. or just Tersoff-Hamman 
    4.5 ! phi in eV                                          
    2 ! Number of voltages to calculate
    1.0 -1.0 ! tip-substrate voltage in V (tip to mass)
    50 ! sampling points in z (perpendicular to the surface)
    10 ! Zmax, maximum tip-surface distance (from 'Zsurf') in \AA
    0.41111 ! Zsurf, origin of the surface in direct coordinates
    0.95555 ! Ztip, origin of the tip in direct coordinates
    0.65000 ! z_s in direct coordinates (neg. values indicate distance from the topmost atom in angs)
    T ! calculate dIdV curve
    -1 1 ! emin, emax for dIdV
    100 ! number of divisions between emin and emax
    2 ! number of points to plot the dIdV curve
    1.721 0.994 16.000 ! coordinates of the point in \AA
    0.000 0.000 16.000 ! coordinates of the point in \AA
    POSCAR ! POSCAR file
    WAVECAR  ! wf file 
    MappingsCAR  !  reciprocal vector and index file
    LGAMMA ! .true. or .false. wether the k-point sampling contains the Gamma point or not.
    T ! do you want the output in cube format?

-----------------------------------------------------------------
                End of input.STMpw
-----------------------------------------------------------------

 Fourth execute STMpw.out

 Output: Bardeen.siesta (unformatted)

 run WSxM and read in Bardeen.siesta

 Last edit
 
    10/09/08  (yr/month/day)
    (2014-07-16) RRR: modification to use vasp4 or vasp5 format.
    (2019-05-24) RRR: modification to use just the TH part.
    (2019-05-28) RRR: modification to support cube files in the TH part.
    (2019-07-16) RRR: modification to include several voltages in the same run.
    (2019-07-17) RRR: modification to calculate dIdV curves.
    (2019-08-07) RRR: modification to calculate MappingsCAR from OUTCAR.
