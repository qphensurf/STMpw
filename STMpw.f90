program STMpw
!
! Simple implementation for Bardeen's transfer hamiltonian
! and tunneling extension for the STM
!
  use declaration
  use determinequantities
  use volumen
  use Fourier
  use currentBardeen
  use mappingscar_gen
  implicit none


! test if input.STMpw exists, if not stop
     inquire (file = 'input.STMpw', exist = fichero)
     if (fichero) then
! reading input.STMpw
     open (8, file = 'input.STMpw', status='old')
     read (8,*, IOSTAT=ios, ERR=100,END=200) phi
     phi = phi / hartree ! all calculations in atomic units
     read (8,*, IOSTAT=ios, ERR=100,END=200) nV
     allocate (V(nV))
     read (8,*, IOSTAT=ios, ERR=100,END=200) (V(jv),jv=1,nV)
     write(*,'(A,30F7.3)') "Voltages to calculate (in V): ",(V(jv),jv=1,nV)
     do jv=1,nV
      V(jv) = V(jv) / hartree ! all calculations in atomic units
     enddo 
     read (8,*, IOSTAT=ios, ERR=100,END=200) N_sampling_z
     read (8,*, IOSTAT=ios, ERR=100,END=200) Zmax
     Zmax = Zmax / bohr
     read (8,*, IOSTAT=ios, ERR=100,END=200) Zsurf
     read (8,*, IOSTAT=ios, ERR=100,END=200) z_s
     read (8,*, IOSTAT=ios, ERR=100,END=200) Bardeen
     if (Bardeen) then
       write(*,'(A)') "Bardeen calculation. NOT TESTED FOR SEVERAL VOLTAGES OR dIdV CURVES!!!"
       read (8,*, IOSTAT=ios, ERR=100,END=200) Ztip
     else
       write(*,'(A)') "Tersoff-Hamman calculation"
       Ztip = 100d0
     endif  
!! matching distance for the tip is hardwired
!! z_t = Ztip - 2.0 Angstroems

! dIdV curve
     read (8,*, IOSTAT=ios, ERR=100,END=200) LDIDV !.true. or .false.
     if (LDIDV) then
       read (8,*, IOSTAT=ios, ERR=100,END=200) Vmin, Vmax
       read (8,*, IOSTAT=ios, ERR=100,END=200) ndiv
       read (8,*, IOSTAT=ios, ERR=100,END=200) npts
       allocate(cdIdV(npts,3))
       do i=1,npts
        read (8,*, IOSTAT=ios, ERR=100,END=200) (cdIdV(i,j),j=1,3)
       enddo 
       if(npts.gt.1000) then
         write(*,'(A)') "We can not calculate more than 1000 dIdV curves."
         npts = 1000
       endif
       cdIdV = cdIdV / bohr
       allocate(V1(ndiv+1))
       allocate(IV(npts,ndiv+1))
       allocate(ngp(npts,3))
       write(*,'(A,F7.3,A,F7.3,A)') "Calculation of IV curves between ",Vmin," eV and ",Vmax," eV"
     endif
     read (8,*, IOSTAT=ios, ERR=100,END=200) NamePOSCAR
     read (8,*, IOSTAT=ios, ERR=100,END=200) NameWF
     read (8,*, IOSTAT=ios, ERR=100,END=200) MAPfile !.true. or .false.
     if (MAPfile) then
      read (8,*, IOSTAT=ios, ERR=100,END=200) NameMAP
      write(*,'(3A)') "We will use '",trim(adjustl(NameMAP)),"' as MappingsCAR file"       
     else
      read (8,*, IOSTAT=ios, ERR=100,END=200) NameOUTCAR
      write(*,'(A)') "We will generate MappingsCAR data from OUTCAR" 
     endif
     read (8,*, IOSTAT=ios, ERR=100,END=200) LGAMMA !.true. or .false.
     read (8,*, IOSTAT=ios, ERR=100,END=200) wsxm   !.true. or .false.
     if(wsxm) read (8,*, IOSTAT=ios, ERR=100,END=200) factor
     read (8,*, IOSTAT=ios, ERR=100,END=200) dat    !.true. or .false.
     read (8,*, IOSTAT=ios, ERR=100,END=200) cube   !.true. or .false.
 
     close (8)
     else
       write(*,*) "You must supply an input.STMpw file."
       stop
     endif

! end of reading input.STMpw

     if(.not.wsxm.AND..not.dat.AND..not.cube.AND..not.LDIDV) then
       write(*,'(A)') "No output has been requested. I refuse to work for nothing."
       stop
     endif  

             call test_POSCAR  !test that matching
! distance is far (<2. Ang) the topmost atom

!determine quantities of the system
        call determine_quantities( NameWF, Number_of_SPIN, Total_Bands , &
     &           Number_of_KPOINTS)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read WF for constant current computations and elastic contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! first reopen with assumed (wrong) record length ICMPLX
      OPEN(UNIT=13,FILE= NameWF,ACCESS='DIRECT', &
                         FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=ICMPLX)

! the first record contains the record length, get it ...
      RDUM=0._q
      READ(13,REC=1,ERR=17421) RDUM,RISPIN ; IDUM=NINT(RDUM)
      IF ((IDUM<=0).OR.(IDUM>10000000)) IDUM=ICMPLX  ! -> error reading
      GOTO 17422
17421 CONTINUE
      IDUM=ICMPLX
17422 CONTINUE
      close (13)
       ISPIN = NINT(RISPIN)
! Allocates 
        ALLOCATE(VKPT(3,Number_of_KPOINTS))
        ALLOCATE (FERTOT (Total_Bands, Number_of_KPOINTS, ISPIN))
        if (MAPfile) allocate (onpl(ISPIN, Number_of_KPOINTS))

! reopen with correct record length (clumsy all that, I know ...)


! Reciprocal vector and wave index file
      unitMAP = 8
! output current units
      unitI = 9
      unitIdat = 19
! output conductance unit
      unitdI = 10
! output conductance unit
      unitdIdat = 25
! output current units Tersoff-Hamman
      unitTH = 11
      unitTHdat = 21
! output conductance unit TH
      unitdITHdat = 20
! output current units Tersoff-Hamman for the TIP
      unitTIP =12
      unitTIPdat =22
! WF input file unit
      unitWF = 13
! output current unit Tersoff-Hamman in cube format
      unitTHcub = 14
! output current unit Tersoff-Hamman in cube format
      unitdITHcub = 15


      OPEN(UNIT=unitWF,FILE=NameWF,ACCESS='DIRECT', &
                      FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IDUM)
! Record 2 read in lattice
      READ(unitWF,REC=2,ERR=200) RKPTSF,RBANDS,ENMAXW, &
     &                         ((A(I,J),I=1,3),J=1,3)
!      write(*,*) "RKPTSF,RBANDS,ENMAXW",RKPTSF,RBANDS,ENMAXW
! transformation into atomic units
! of the unit cell vectors
           A = A /bohr
! generate reciprocal unit cell

! third dimension in 90 degrees with respect to surface
        if (A (1,3) > 0.1 .or. A(2,3) > 0.1 ) then

              print *, 'The calculation has to be repeated'
              print *, 'with the 3rd axis perpendicular to the surface plane!'
              stop

        endif
! test on cells to check that ionic positions are the same
! for WAVECAR and POSCAR 
! we test that the 3rd axis is the same withing 0.1 Angstroems
!
        if (abs (cell(3,3)*lat_par - A(3,3) * bohr) > 0.1 ) then
            print *, 'The cells in POSCAR and WAVECAR'
            print *, 'are more than 0.1 angstroems different!!'
            print *, 'STOP and check the input files!'
            stop
        endif
       
      call volume (A(:,1),A(:,2),A(:,3),vol)
      call surface (A(:,1),A(:,2),surf)
!     call reciprocal (A(:,1),A(:,2),A(:,3),vol,vec(:,1),vec(:,2),vec(:,3))

! When reading it from MappingsCAR, you don't need to calculate the
! reciprocal vector:
! Reading MappingsCAR
! Reciprocal vector and wave index file
   
      if(.not.MAPfile) then
       call mappings_from_outcar(NameOUTCAR,ENMAXW, ONPL, oerror)
       if(oerror.gt.0) then
         write(*,'(A)') "Generation of OUTCAR failed"
         stop
       endif  
       NameMAP = "MappingsCAR_gen"  
      endif

      open (unitMAP, file = NameMAP, form = 'unformatted' )      
      read (unitMAP) Efermi 
      write(*,'(A,F9.4,A)') "Fermi energy at ",Efermi," eV"
           Efermi = Efermi / hartree ! all calculations in atomic units
      read (unitMAP) NSPIN, NKPTS, NGX, NGY, NGZ

! general allocation
       call allocation

! spin factor 
      if (Number_of_SPIN == 1) then
             spin_factor = 2
      else
             spin_factor = 1
      endif


      if (NSPIN /= Number_of_SPIN) then
          print *, 'WAVECAR and MappingsCAR do not match: '
          print *, ' different spin.'
          stop
      endif
      if (NKPTS /= Number_of_KPOINTS) then
          print *, 'WAVECAR and MappingsCAR do not match: '
          print *, ' different number of K points.'
          stop
      endif


      Zmin = (z_s - Zsurf) * A(3,3)

      stepX = sqrt(sum(A(:,1)**2))/(Ngx)
      stepY = sqrt(sum(A(:,2)**2))/(Ngy)
! matching point to tip distance
      distance = (Ztip-z_t) * A(3,3)

! z is the difference
! between matching points in this program
! in other words: when z=0 the distance
! between tip and sample is the offset
      offset = Zmin + distance
      stepZ = Zmax / (N_sampling_z-1)
! Hence the initial distance is "offset"
! the maximum distance is "Zmax+offset"

      if (LDIDV) then
!       do i=1,npts
!         if(cdIdV(i,3).lt.z_s*A(3,3)) then
!           write(*,'(A)') "We do not calculate the dIdV curves because one point is out of range:"
!           write(*,*) cdIdV(i,3)*bohr," < ", z_s*A(3,3)*bohr
!           LDIDV = 0
!         elseif(cdIdV(i,3).gt.(z_s*A(3,3)+Zmax)) then
!           write(*,'(A)') "We do not calculate the dIdV curves because one point is out of range:"
!           write(*,*) cdIdV(i,3)*bohr," > ", (z_s*A(3,3)+Zmax)*bohr
!           LDIDV = 0
!         endif
!       enddo

        do jv=1,ndiv+1
          V1(jv) = Vmin + (Vmax-Vmin)/(ndiv-1)*(jv-1) 
          V1(jv) = V1(jv) / hartree
!          write(*,*) "V(",jv,")=",V1(jv)
        enddo
        ndiv = ndiv + 1

        do iy = 0, ngy
          do ix = 0, ngx
            cngx(ix,iy)=cell(1,1)*ix/(ngx+1)+cell(2,1)*iy/(ngy+1)
            cngy(ix,iy)=cell(1,2)*ix/(ngx+1)+cell(2,2)*iy/(ngy+1)
            cngx(ix,iy)=cngx(ix,iy)*lat_par
            cngy(ix,iy)=cngy(ix,iy)*lat_par
          enddo
        enddo
!
        write(*,'(A)') "Points to calculate dIdV curves (in angs):"
        do i=1,npts
         ngp(i,1)=ngx+100
         do iy = 0, ngy
         do ix = 0, ngx
          if(abs(cdIdV(i,1)-cngx(ix,iy)/bohr).lt.stepX/2.AND.abs(cdIdV(i,2)-cngy(ix,iy)/bohr).lt.stepY/2) then
           ngp(i,1)=ix
           ngp(i,2)=iy
          endif
         enddo
         enddo
          if(ngp(i,1).eq.ngx+100) then
           write(*,'(A,i4,A)') "We do not calculate the dIdV curves because point",i," is out of range:"
           write(*,'(A,2F8.3,A)') "(",(cdIdV(i,j)*bohr,j=1,2),") not in the unit cell."
           LDIDV = .false.
          endif
         ngp(i,3)=N_sampling_z+100
         do iz=1,N_sampling_z
          if(abs(cdIdV(i,3)-(z_s*A(3,3)+stepZ*(iz-1))).lt.stepZ/2) then
           ngp(i,3)=iz
          endif
         enddo
          if(ngp(i,3).eq.N_sampling_z+100) then
           write(*,'(A,i4,A)') "We do not calculate the dIdV curves because point",i," is out of range:"
           write(*,'(F8.3,A,F8.3,A,F8.3,A)') & 
                   cdIdV(i,3)*bohr," not in [", z_s*A(3,3)*bohr,",",(z_s*A(3,3)+Zmax)*bohr,"]"
           LDIDV = .false.
          endif
        write(*,'(A,i3,A,3F8.3,A)') " Requested point",i," = (",(cdIdV(i,j)*bohr,j=1,3),")"
        write(*,'(A,i3,A,3F8.3,A)') "    Actual point",i," = (",&
           cngx(ngp(i,1),ngp(i,2)),cngy(ngp(i,1),ngp(i,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(i,3)-1)),")"
        enddo

      endif  

!Initialise records
      IREC = 2

! compute current as we read WF file

     if(Bardeen) then
       intensity = 0
       dintensity_dV = 0
       if(LDIDV) intensity2 = 0
       Tersoff_t =  0
     endif
     Tersoff_c =  0
     Tersoff_s =  0
     if(LDIDV) Tersoff_s2 = 0

!     density_z = 0

!Loop on spin and k-points

     SPIN:  do spin1=1, Number_of_SPIN
      
      write(*,'(A,i2,A)') " Calculating spin ", spin1, " ..."

         KPOINTS: do kp=1, Number_of_KPOINTS

          write(*,'(A,i4,A)') "  Calculating kpoint ", kp, " ..."

!read in orbitals

         IREC = IREC + 1

!  eigen energies:

       read (unitWF, REC = IREC) RNPL,VKPT (1,kp), &
                        VKPT (2,kp), VKPT (3,kp), &
                       (EIG(J),FERTOT(J,kp,spin1),J=1,Total_Bands)

           EIG = EIG/hartree - Efermi

        NPL=NINT(RNPL)

        if(.not.MAPfile) then
          if(NPL.ne.onpl(spin1,kp)) then
           write(*,'(A)') "The number of planewaves is not the same in WAVECAR and MappingsCAR"
           write(*,'(A)') " We can not continue. It might be a problem reading ENCUT from OUTCAR."
           write(*,'(A,I9,A,I9)') " NPL from WAVECAR = ", NPL, ", ONPL from MappingsCAR = ", onpl(spin1,kp)
           stop
          endif
        endif 
        
       allocate (IGX(NPL),IGY(NPL),IGZ(NPL))
! reciprocal vectors
       i = 0
       do iy = 1, NGY
        do ix = 1, NGX
         i = i+1
         read (unitMAP) gx (i), gy (i)
        enddo
       enddo
! reading wave function indices
       do i = 1, NPL
        read (unitMAP) IGX (i), IGY (i), IGZ (i)
!              IGX (i) = mod (IGX (i) + ngx/2, ngx)
!              IGY (i) = mod (IGY (i) + ngy/2, ngy)
       enddo

! change to atomic units
     gx = gx * bohr
     gy = gy * bohr

     write(*,*) "   Reading and matching wavefunctions ..."

     t_wf = 0.0_q
     t_ma = 0.0_q

     call cpu_time(start)

         BANDS: do iband=1,Total_Bands

       IREC = IREC + 1

       CW =(0.0,0.0)
       call cpu_time(st_tmp)
       read (unitWF, REC = IREC) ( CW (I), I = 1, NPL )
       call cpu_time(fi_tmp)
       t_wf = t_wf + fi_tmp - st_tmp

       
!uncomment test on the electron density to check
! that everything is OK
! density on the spot x =0 y =0 and dependence on z
!BEGIN test
!       call density(density1,CW,ngx, ngy, ngz, IGX,IGY,IGZ,NPL)
!       density_z(:) = density_z(:) + Fermi(dreal(EIG(iband)))*density1(:)/Number_of_KPOINTS
!MID test


! obtain 2-D coefficients for the surface and the tip
! by performing the exponential matching

   call cpu_time(st_tmp)
   call matching (CW, A_S, C_T, ngx, ngy, ngz, z_s,z_t,IGX,IGY,IGZ,NPL,temp)
   call cpu_time(fi_tmp)
   t_ma = t_ma + fi_tmp - st_tmp

     A_G(:,:,iband)=A_S(:,:)
     C_G(:,:,iband)=C_T(:,:)
! close Loop on bands

           enddo BANDS
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "     Time to read and match wavefunctions = ", finish-start, " s"
     write (*,'(A, F8.1, A)') "       Time to read wavefunctions = ", t_wf, " s"
     write (*,'(A, F8.1, A)') "       Time to match wavefunctions = ", t_ma, " s"
!END test
!     do iz = 1, ngz
!     print *, iz, density_z(iz)/vol/(bohr*bohr*bohr)
!     enddo
!     stop

! Compute the current calculation for all different tip
! positions in the cell and from Z_min to Z_max.
! Loop on tip-surface distances:

       call cpu_time(start)
       write(*,*) "   Computing current ..."

! Loop on voltages
   
   VOLTAGES: do jv = 1, nV
      
       write(*,'(A,F7.3)') "      Computing V = ", V(jv)*hartree

       Heights: do iz=1,N_sampling_z

       z=stepZ*(iz-1)


       if (Bardeen) then
! Loop on tip bands
         do tband=1,Total_Bands


! Loop on substrate bands
         do sband=1,Total_Bands
  
! Check voltage sign (tip to mass, then positive
! means empty substrate states)
   
    if ( V(jv) < 0 ) then

          Fermi_t = 1-Fermi(dreal(EIG(tband)))
          Fermi_s = Fermi (dreal(EIG(sband)))

     constant_I= - 2 * pi_d /( sqrpi * sigma)
     constant_dIdV=  4 * pi_d /( sqrpi * sigma**3)

    else

          Fermi_t = Fermi(dreal(EIG(tband)))
          Fermi_s = 1-Fermi (dreal(EIG(sband)))

     constant_I=  2 * pi_d /( sqrpi * sigma)
     constant_dIdV= - 4 * pi_d /( sqrpi * sigma**3)

    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     constant_dIdV=spin_factor*constant_dIdV*surf**2/vol**2

! before calculating we check for empty states ...
! and energy conservation

             if (Fermi_t*Fermi_s > 0.1) then

             delta_E = ((dreal(EIG(tband)-EIG(sband)) + V(jv))/sigma) **2

             if (delta_E <3.) then

! compute Bardeen's matrix element: currentSQ which is already
! the squared modulus

     A_S(:,:)=A_G(:,:,sband)
     C_T(:,:)=C_G(:,:,tband)
 
   call matrix_element (A_S,C_T,ngx,ngy,phi,gx,gy,z,currentSQ)

! compute current

      if (.not.LGAMMA) then
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (Number_of_KPOINTS)
      elseif (kp == 1) then
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      else
   intensity(:,:,iz,jv) = intensity (:,:,iz,jv) + 2.0*constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      end if



! compute conductance

      if (.not.LGAMMA) then
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (2*Number_of_KPOINTS-1)
      else
dintensity_dV(:,:,iz,jv) = dintensity_dV(:,:,iz,jv) + 2.0*constant_dIdV * Fermi_t * &
    & Fermi_s  * currentSQ(:,:) * exp (-delta_E) *  &
    & (dreal(EIG(tband)-EIG(sband)) + V(jv))/ &
       & (2*Number_of_KPOINTS-1)
      end if

             endif
             endif

          enddo
          enddo

      endif ! Bardeen part

! Redo calculation to
! calculate Tersoff-Hamman picture

    if ( V(jv) < 0 ) then
     constant_I= - 2 * pi_d /( sqrpi * sigma)
    else
     constant_I=  2 * pi_d /( sqrpi * sigma)
    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     
! Loop on bands
         do tband=1,Total_Bands

         if ( V(jv) > 0) then
         if (dreal(EIG(tband)) > 0 .and. dreal(EIG(tband)) < V(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (Number_of_KPOINTS)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I 
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + 2.0*TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

        endif

        else

        if (dreal(EIG(tband)) < 0 .and. dreal(EIG(tband)) > V(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (Number_of_KPOINTS)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_c(:,:,iz,jv) =  Tersoff_c(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1) * exp (-((dreal(EIG(tband)-V(jv))/sigma)**2))&
       & * constant_I
         if(Bardeen) Tersoff_t(:,:,iz,jv) =  Tersoff_t(:,:,iz,jv) + 2.0*TH_t (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
         Tersoff_s(:,:,iz,jv) =  Tersoff_s(:,:,iz,jv) + 2.0*TH_s (:,:)/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

          endif
          endif

! Loop on bands
         enddo

 
           enddo Heights
         enddo VOLTAGES  
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "    Time for current = ", finish-start, " s"

! Calculation on dIdV curves

    if(LDIDV) then 
       call cpu_time(start)
       write(*,*) "   Computing dIdV curves ..."

  npts_dIdV: do ip = 1,npts
   dIdV: do jv = 1, ndiv
      
       iz=ngp(ip,3)

       z=stepZ*(iz-1)


       if (Bardeen) then
! Loop on tip bands
         do tband=1,Total_Bands


! Loop on substrate bands
         do sband=1,Total_Bands
  
! Check voltage sign (tip to mass, then positive
! means empty substrate states)
   
    if ( V1(jv) < 0 ) then

          Fermi_t = 1-Fermi(dreal(EIG(tband)))
          Fermi_s = Fermi (dreal(EIG(sband)))

     constant_I= - 2 * pi_d /( sqrpi * sigma)
     constant_dIdV=  4 * pi_d /( sqrpi * sigma**3)

    else

          Fermi_t = Fermi(dreal(EIG(tband)))
          Fermi_s = 1-Fermi (dreal(EIG(sband)))

     constant_I=  2 * pi_d /( sqrpi * sigma)
     constant_dIdV= - 4 * pi_d /( sqrpi * sigma**3)

    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     constant_dIdV=spin_factor*constant_dIdV*surf**2/vol**2

! before calculating we check for empty states ...
! and energy conservation

             if (Fermi_t*Fermi_s > 0.1) then

             delta_E = ((dreal(EIG(tband)-EIG(sband)) + V1(jv))/sigma) **2

             if (delta_E <3.) then

! compute Bardeen's matrix element: currentSQ which is already
! the squared modulus

     A_S(:,:)=A_G(:,:,sband)
     C_T(:,:)=C_G(:,:,tband)
 
   call matrix_element (A_S,C_T,ngx,ngy,phi,gx,gy,z,currentSQ)

! compute current

      if (.not.LGAMMA) then
   intensity2(ip,jv) = intensity2 (ip,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (Number_of_KPOINTS)
      elseif (kp == 1) then
   intensity2(ip,jv) = intensity2 (ip,jv) + constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      else
   intensity2(ip,jv) = intensity2 (ip,jv) + 2.0*constant_I * Fermi_t *   &
    & Fermi_s  * currentSQ(ngp(ip,1),ngp(ip,2)) * exp (-delta_E)/                    &
    & (2*Number_of_KPOINTS-1)
      end if

             endif
             endif

          enddo
          enddo

      endif ! Bardeen part

! Redo calculation to
! calculate Tersoff-Hamman picture

    if ( V1(jv) < 0 ) then
     constant_I= - 2 * pi_d /( sqrpi * sigma)
    else
     constant_I=  2 * pi_d /( sqrpi * sigma)
    endif
     constant_I=spin_factor*constant_I*surf**2/vol**2
     
! Loop on bands
         do tband=1,Total_Bands

         if ( V1(jv) > 0) then
         if (dreal(EIG(tband)) > 0 .and. dreal(EIG(tband)) < V1(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + 2.0*TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if


         endif

         else

         if (dreal(EIG(tband)) < 0 .and. dreal(EIG(tband)) > V1(jv) ) then 
      
     A_S(:,:)=A_G(:,:,tband)
     C_T(:,:)=C_G(:,:,tband)

   call matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z,TH_t,TH_s)

      if (.not.LGAMMA) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (Number_of_KPOINTS)
      elseif (kp == 1) then
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      else
         Tersoff_s2(ip,jv) =  Tersoff_s2(ip,jv) + 2.0*TH_s (ngp(ip,1),ngp(ip,2))/vol/ &
       & (2*Number_of_KPOINTS-1)
      end if

          endif
          endif

! Loop on bands
         enddo

         enddo dIdV      
        enddo npts_dIdV
     call cpu_time(finish)
     write (*,'(A, F8.1, A)') "    Time for computing dIdV curves = ", finish-start, " s"
    endif 

      deallocate (IGX,IGY,IGZ)

         enddo KPOINTS
       enddo SPIN

      if(.not.MAPfile) then
       cmd = 'rm -f MappingsCAR_gen'
       call system(cmd)
      endif 


! output in WSxM format
! *.siesta
! output in ASCII to be converted into STM.dato
! and read by imagen.f into Mathematica

        write(*,*) "Writing files ..."
        call cpu_time(start)
        if(wsxm) write(*,'(A,F8.0)') "Output in WSxM format. Output will be multiplied by ",factor
        if(dat) write(*,'(A)') "Output for gnuplot"
        if(cube) write(*,'(A)') "Output for TH in cube format"
      do jv=1,nV 
        write(volt,'(F6.3)') V(jv)*hartree
        cmd = 'mkdir V_'//trim(adjustl(volt))
        call system(cmd) 
        UCELL = A
        UCELL (3,3) = Zmax

        if(wsxm) then
         name_file = 'V_'//trim(adjustl(volt))//'/Bardeen_V_'//trim(adjustl(volt))//'.siesta'
         if (Bardeen) open (unitI,file=name_file, form = 'unformatted')
         name_file = 'V_'//trim(adjustl(volt))//'/TH_V_'//trim(adjustl(volt))//'.siesta'
         open (unitTH,file=name_file, form = 'unformatted')
         name_file = 'V_'//trim(adjustl(volt))//'/TH_tip_V_'//trim(adjustl(volt))//'.siesta'
         if (Bardeen) open (unitTIP,file=name_file, form = 'unformatted')
         if (Bardeen) write (unitI) UCELL
         write (unitTH) UCELL
         if (Bardeen) write (unitTIP) UCELL
         if (Bardeen) write (unitI) ngx, ngy, N_sampling_z,1
         write (unitTH) ngx, ngy, N_sampling_z,1
         if (Bardeen) write (unitTIP) ngx, ngy, N_sampling_z,1
         do iz= 1, N_sampling_z
          do iy = 0, ngy-1
           if (Bardeen) write (unitI) (real(abs(intensity (ix,iy,iz,jv))*factor), ix=0,ngx-1)
           write (unitTH) (real(Tersoff_s(ix,iy,iz,jv)*factor), ix=0,ngx-1)
           if (Bardeen) write (unitTIP) (real(Tersoff_t(ix,iy,iz,jv)*factor), ix=0,ngx-1)
          enddo
         enddo
         if (Bardeen) close (unitI)
         close (unitTH)
         if (Bardeen) close (unitTIP)

! Output for conductance.
! Output in WSxM format,
! not extremely useful because you can't do a scan on voltage
! but you can get a conductance map

!         if (Bardeen) then
!           name_file = 'V_'//trim(adjustl(volt))//'/dIdV_Bardeen_V_'//trim(adjustl(volt))//'.siesta'
!           open (unitdI,file=name_file, form = 'unformatted')
!           write (unitdI) UCELL
!           write (unitdI) ngx, ngy, N_sampling_z,1
!           do iz= 1, N_sampling_z
!             do iy= 0, ngy-1
!               write (unitdI) (real(dintensity_dV (ix,iy,iz,jv)*factor), ix=0,ngx-1)
!             enddo
!           enddo
!           close (unitdI)
!         endif
        endif

        if(dat) then
          name_file = 'V_'//trim(adjustl(volt))//'/Bardeen.dat'
          if (Bardeen) open (unitIdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/TH.dat'
          open (unitTHdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/TH_tip.dat'
          if (Bardeen) open (unitTIPdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_TH.dat'
          open (unitdITHdat,file=name_file)
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_Bardeen.dat'
          if (Bardeen) open (unitdIdat,file=name_file)
          write(unitdITHdat,*) N_sampling_z, ngy, ngx
          write(unitTHdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitTIPdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitIdat,*) N_sampling_z, ngy, ngx
          if (Bardeen) write(unitdIdat,*) N_sampling_z, ngy, ngx
          do iz= 1, N_sampling_z
           do iy = 0, ngy-1
            do ix = 0, ngx-1
             write (unitdITHdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_c(ix,iy,iz,jv)
             write (unitTHdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_s(ix,iy,iz,jv)
             if (Bardeen) write (unitTIPdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)), bohr*(iy)*stepY, bohr*(ix)*stepX, Tersoff_t(ix,iy,iz,jv) 
             if (Bardeen) write (unitIdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, intensity (ix,iy,iz,jv)
             if (Bardeen) write (unitdIdat,'(4g14.4)') &
               bohr*(Zmin+stepZ*(iz-1)),bohr*(iy)*stepY, bohr*(ix)*stepX, dintensity_dV (ix,iy,iz,jv)
            enddo
           enddo
          enddo

          close(unitTHdat)
          if (Bardeen) close(unitTIPdat)
          if (Bardeen) close(unitIdat)
          if (Bardeen) close(unitdIdat)
          close(unitdITHdat)
        endif

        if (cube) then
          name_file = 'V_'//trim(adjustl(volt))//'/TH_V_'//trim(adjustl(volt))//'.cube'
          open (unitTHcub,file=name_file)
          write(unitTHcub,'(A)') "STM image in TH approximation."
          write(unitTHcub,'(3A,F6.3,A)') "V = ",volt," V. z = 0 is at",Zmin*bohr," angs over the surface."
          write(unitTHcub,'(i5,3F12.6)') nat, 0.0, 0.0, z_s*cell(3,3)*lat_par/bohr
          write(unitTHcub,'(i5,3F12.6)') ngx, (UCELL(j,1)/ngx,j=1,3)
          write(unitTHcub,'(i5,3F12.6)') ngy, (UCELL(j,2)/ngy,j=1,3)
          write(unitTHcub,'(i5,3F12.6)') N_sampling_z, (UCELL(j,3)/N_sampling_z,j=1,3)
          do i=1,nat
            write(unitTHcub,'(i5,4F12.6)') zat(i),0.000,(coord(i,j)/bohr,j=1,3)
          enddo  
          do ix = 0, ngx-1
            do iy = 0, ngy-1
              write(unitTHcub, '(6E13.5)') (Tersoff_s(ix,iy,iz,jv), iz=1,N_sampling_z)
            enddo
          enddo  
          close(unitTHcub)
! dIdV          
          name_file = 'V_'//trim(adjustl(volt))//'/dIdV_TH_V_'//trim(adjustl(volt))//'.cube'
          open (unitdITHcub,file=name_file)
          write(unitdITHcub,'(A)') "dIdV in TH approximation."
          write(unitdITHcub,'(3A,F6.3,A)') "V = ",volt," V. z = 0 is at",Zmin*bohr," angs over the surface."
          write(unitdITHcub,'(i5,3F12.6)') nat, 0.0, 0.0, z_s*cell(3,3)*lat_par/bohr
          write(unitdITHcub,'(i5,3F12.6)') ngx, (UCELL(j,1)/ngx,j=1,3)
          write(unitdITHcub,'(i5,3F12.6)') ngy, (UCELL(j,2)/ngy,j=1,3)
          write(unitdITHcub,'(i5,3F12.6)') N_sampling_z, (UCELL(j,3)/N_sampling_z,j=1,3)
          do i=1,nat
            write(unitdITHcub,'(i5,4F12.6)') zat(i),0.000,(coord(i,j)/bohr,j=1,3)
          enddo  
          do ix = 0, ngx-1
            do iy = 0, ngy-1
              write(unitdITHcub, '(6E13.5)') (abs(Tersoff_c(ix,iy,iz,jv)), iz=1,N_sampling_z)
            enddo
          enddo  
          close(unitdITHcub)
        endif 

     enddo ! Voltages

     if (LDIDV) then
      cmd = 'mkdir dIdV_curves'
      call system(cmd)
      do ip=1,npts
        unitIV=1000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/IV_TH_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# IV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "# Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         if(V1(jv).lt.0) Tersoff_s2(ip,jv) = -1.0 * Tersoff_s2(ip,jv)
         write(unitIV,*) V1(jv)*hartree,Tersoff_s2(ip,jv)
        enddo
        close (unitIV)
      enddo ! points for dIdV
      do ip=1,npts
        unitIV=2000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/dIdV_TH_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# dIdV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "#   Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         write(unitIV,*) V1(jv)*hartree,(Tersoff_s2(ip,jv+1)-Tersoff_s2(ip,jv))/(V1(jv+1)-V1(jv))/hartree
        enddo
        close (unitIV)
      enddo ! points for dIdV
!
     if (Bardeen) then
      do ip=1,npts
        unitIV=3000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/IV_Bardeen_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# IV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "# Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         if(V1(jv).lt.0) intensity2(ip,jv) = -1.0 * intensity2(ip,jv)
         write(unitIV,*) V1(jv)*hartree,intensity2(ip,jv)
        enddo
        close (unitIV)
      enddo ! points for dIdV
      do ip=1,npts
        unitIV=4000+ip
        write(cnpt,'(I6)') ip
        name_file = 'dIdV_curves/dIdV_Bardeen_npt_'//trim(adjustl(cnpt))//'.dat'
        open (unitIV,file=name_file)
!        write(unitIV,'(A,3F7.3,A)') "# dIdV for point (",(cdIdV(ip,j)*bohr,j=1,3),")"
        write(unitIV,'(A,3F7.3,A)') "#   Actual point (", &
                cngx(ngp(ip,1),ngp(ip,2)),cngy(ngp(ip,1),ngp(ip,2)),bohr*(z_s*A(3,3)+stepZ*(ngp(ip,3)-1)),")"
        do jv=1,ndiv-1
         write(unitIV,*) V1(jv)*hartree,(intensity2(ip,jv+1)-intensity2(ip,jv))/(V1(jv+1)-V1(jv))/hartree
        enddo
        close (unitIV)
      enddo ! points for dIdV
     endif ! Bardeen
    endif ! dIdV curves

     call cpu_time(finish)
     write (*,'(A, F8.1, A)') " Time to write files = ", finish-start, " s"

100  if (ios>0) then
     print *,'ERROR reading input.STMpw (wrong format or corrupt file), IOS=',ios
     endif
200  if (ios<0) then
     print *,'End of input.STMpw (empty file or not enough data), IOS=',ios
     endif

CONTAINS

! allocation of memory

       subroutine allocation

       allocate (EIG(Total_bands))
       allocate (A_G(0:ngx-1,0:ngy-1,Total_bands))
       allocate (C_G(0:ngx-1,0:ngy-1,Total_bands))
       allocate (A_S(0:ngx-1,0:ngy-1))
       allocate (C_T(0:ngx-1,0:ngy-1))
       allocate (cngx(0:ngx,0:ngy))
       allocate (cngy(0:ngx,0:ngy))
       allocate (CW(ngx*ngy*ngz))
!       allocate (density_z(ngz))
!       allocate (density1(ngz))
       allocate (temp(0:ngx-1,0:ngy-1,0:ngz-1))

       if(Bardeen) then
         allocate (intensity(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         allocate (Tersoff_t(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         allocate (dintensity_dV(0:ngx-1,0:ngy-1,N_sampling_z,nV))
         if(LDIDV) allocate (intensity2(npts,ndiv))
       endif  
       allocate (Tersoff_s(0:ngx-1,0:ngy-1,N_sampling_z,nV))
       allocate (Tersoff_c(0:ngx-1,0:ngy-1,N_sampling_z,nV))
       if(LDIDV) allocate (Tersoff_s2(npts,ndiv))
       allocate (currentSQ(0:ngx-1,0:ngy-1))
       allocate (TH_s(0:ngx-1,0:ngy-1))
       allocate (TH_t(0:ngx-1,0:ngy-1))
       allocate (gx(NGX*NGY))
       allocate (gy(NGX*NGY))

       return
       end subroutine allocation

! test that matching
! distance is far (<2. Ang)eng) the topmost atom
      subroutine test_POSCAR 
      type vasp
      real (q), dimension(3)       :: coord
      character( len = 1 ), dimension (3)  :: species
      end type vasp

      character( len = 80 ) :: line
      character( len = 1 ) :: charac
      character( len = 10 ) :: din
      
      integer :: n_atom_species(90) = 0
      real (q) :: maximum_surf, maximum_tip, new_coor, coor_max, coor_tip, new_tip, c(3)
!      real (q), allocatable :: coord(:,:)
      integer :: ii, i1, i2

      type( vasp )   vasp_line

      logical :: fractional, selective = .false.

      character (len = 2) :: typelist(50) = "0", norm_name
      integer :: ntyp1, znumber(50) = 0, nionlist(80) = 0
!      integer, allocatable :: zat(:)

      integer, parameter  :: nel = 103
      character(len=2), parameter, dimension(nel) :: atname = &
               (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
                 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
                 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
                 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                 'Md','No','Lr'/)

      maximum_surf = 0.0
      maximum_tip = 0.0

      open (8, file = NamePOSCAR)


      read(8,*)
      read(8,*) lat_par
      read(8,*) cell(1,1), cell(1,2), cell(1,3)
      read(8,*) cell(2,1), cell(2,2), cell(2,3)
      read(8,*) cell(3,1), cell(3,2), cell(3,3)

        read (8,'(A)') line

        read(line,*,end=999) charac

        if(charac>='0' .AND. charac<='9') then
          write(*,'(A)') 'Only Vasp 5 format is supported'
          stop
        else
          write(*,'(A)') 'Vasp 5 format'
        endif
!
      i=1
      j=1
      do
        if (i.eq.80) exit
        if (line(i:i).ne." ") then
          typelist(j) = norm_name(line(i:i+1))
!          write(*,"(I4,2A)") j," ",typelist(j)
          i = i+2
          j = j+1
        else
          i = i+1
        endif
      enddo
      ntyp1 = j-1
!
      do i = 1,ntyp1
        do j = 1,nel
          if (atname(j) .eq. typelist(i)) then
            znumber(i) = j
!            write(*,*) typelist(i),znumber(i)
            exit
          endif
        enddo
        if(znumber(i).eq.0) then
          write(*,'(A)') "Can not identify element ",typelist(i)
          stop
        endif
      enddo

        read(8,'(A)') line
        read(line,*,end=999) n_atom_species
999     continue
        do i=1,90
         if(n_atom_species(i).eq.0) exit
        enddo
        n_species=i-1        
        if (n_species .ne. ntyp1 ) then
          write(*,'(A)') "POSCAR file is inconsistent"
          stop
        endif      
!        write(*,'(i3,A)') n_species," n_species."

      nat = sum(n_atom_species)

      allocate (coord(nat, 3))
      allocate (zat(nat))

      read(8,'(A)') line
      read(line,*,end=999) charac
!
      if(charac.eq.'S' .OR. charac.eq.'s') then
       write(*,'(A)') "Selective dynamics"
       read(8,'(A)') line
      endif
      read(line,*,end=999) charac
!
      if(charac.eq.'C' .OR. charac.eq.'c'.OR. charac.eq.'K'.OR. charac.eq.'k') then
       write(*,'(A)') "We found cartesian coordinates, but only fractional coordinates are supported."
       stop
      else
       write(*,'(A)') "Fractional coordinates"
       fractional = .true.
      endif

!reading coordinates
      ii=0
      do i1=1,ntyp1
       do i2=1,n_atom_species(i1)
        ii=ii+1
        if( selective ) then
          read(8,*) (coord(ii,j),j=1,3),din
        else
          read(8,*) (coord(ii,j),j=1,3)
        endif
        zat(ii)=znumber(i1)
!        write(*,300) (c(j),j=1,3),i1,ii,typelist(i1),znumber(i1)
!        do j=1,3
!          c(j)=cell(1,j)*coord(ii,1)+cell(2,j)*coord(ii,2)+cell(3,j)*coord(ii,3)
!          c(j)=c(j)*lat_par
!        enddo
!        write(*,300) (c(j),j=1,3),i1,ii,typelist(i1),znumber(i1)
       enddo
      enddo
      
300   format(3F20.14,2(I4,1x),A4,2i4)

! test on distances
      do j = 1, nat

!               read(8,*) vasp_line
         
!               new_coor = mod(vasp_line%coord(3)+1,1._q)
               new_coor = mod(coord(j,3)+1,1._q)
               new_tip = new_coor

               if (new_coor < 0.5*(Zsurf+Ztip) .and. new_coor >= Zsurf) then

                 if (new_coor-Zsurf > maximum_surf) then
                           maximum_surf = new_coor-Zsurf
                           coor_max = new_coor
                 endif

               endif
! tip matching distance
                if (new_tip < Ztip .and. new_tip >= 0.5*(Zsurf+Ztip)) then
                  if (Ztip-new_tip > maximum_tip) then
                          maximum_tip = Ztip-new_tip
                          coor_tip = new_tip
                  endif
                endif


      end do

      write(*,'(A,F10.5)') "Coord. max. = ",coor_max

! we define z_t (matching for the tip here:
! z_t = Ztip - 2.0 angstroems

       z_t = coor_tip - 2.0 / (lat_par * cell(3,3))

     if(Bardeen) then
      write(*,'(A,F10.5)') "Coord. tip = ", coor_tip
      write(*,'(A,F10.6)') "z_t = ", z_t
     endif

! if z_s < 0 it indicates the distance above the topmost atom
      if (z_s < 0) then
        z_s = coor_max - z_s / (lat_par * cell(3,3))
        write(*,'(A,F10.6)') "z_s = ", z_s
      endif  

! if using defaults
! calculate z_s as 2.5 Ang above the topmost atom

      if (.not.fichero) then

       z_s = coor_max + 2.5 / (lat_par * cell(3,3))

      end if

      if (Bardeen.AND.(z_s > z_t)) then
         print *, 'WARNING:'
         print *, ' z_s MUST be smaller than z_t!'
         print *, 'WARNING: we do not stop but the STM calculation is wrong!!!'
         print *, 'the input data are:'
         print *, 'Ztip (tip base)=',Ztip
         print *, 'z_t (matching for tip)=',z_t
         print *, 'tip matching distance from the tip base is=', (Ztip-z_t)**lat_par*cell (3,3)
         print *, 'Zsurf (surface base)=',Zsurf
         print *, 'z_s (matching for surface)=',z_s
      endif

            if (coor_max*lat_par*cell (3,3) + 2.0 > z_s *lat_par*cell (3,3)) then
                 print *, 'WARNING: Surface topmost atom (Ang) =', (coor_max-Zsurf)*lat_par*cell (3,3)
                 print *, 'WARNING: z_s-Zsurf (Ang) =', (z_s-Zsurf) *lat_par*cell (3,3)
                 print *, 'WARNING: your matching distance from the surface is less'
                 print *, 'WARNING: than 2 Angstroems from the topmost atom!!!'
         print *, 'tip matching distance from the tip base is (Ang)=', (Ztip-z_t)**lat_par*cell (3,3)
                 print *, 'WARNING: we do not stop but the STM calculation is wrong!!!'
            endif


         close (8)
      if (Zsurf > z_s) then
         print *, ' Zsurf MUST be smaller than z_s!'
         print *, ' No output.'
         stop
      endif
      do i=1,nat
        do j=1,3
          c(j)=cell(1,j)*coord(i,1)+cell(2,j)*coord(i,2)+cell(3,j)*coord(i,3)
          c(j)=c(j)*lat_par
        enddo
        do j=1,3
          coord(i,j)=c(j)  
        enddo
      enddo  
 
       return
       end subroutine test_POSCAR 


end program STMpw
