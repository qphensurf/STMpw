!
! Program to generate MappingsCAR from OUTCAR
!
! Roberto Robles, 2019/08/07
!
      module mappingscar_gen
              
      contains        

      subroutine mappings_from_outcar(NameOUTCAR,enmaxw, onpl, oerror)
              
      implicit none        

      double precision :: enmaxw, enmax, efermi, cell_r(3,3), etotal
      double precision :: g1, g2, g3, gix, giy, giz
      double precision, allocatable :: kpt(:,:)

      integer :: i, j, ix, iy, iz, is, ik, oerror, isym
      integer :: nspin, nkpts, nbands, ngx, ngy, ngz, status
      integer, allocatable :: rangex(:), rangey(:), rangez(:) 
      integer, pointer :: onpl(:,:)

      character( len = 150) :: line, name

      real, parameter :: tpi = 6.28318530717958647693d0
      real, parameter :: c = 0.26246582250210965422d0

      logical :: present, cell = .false., lfermi = .false.
      logical :: lngx = .false.

      character ( len = 99) :: NameOUTCAR

! extract parameters from OUTCAR

      write(*,'(A)') "Generating MappingsCAR from OUTCAR:"

      inquire(file=NameOUTCAR, exist=present)
      if(.NOT.present) then
       write(*,'(2A)') NameOUTCAR,' is not present'
       oerror = 1
       return
      endif

      open(29, file=NameOUTCAR)

      do
       read(29,'(A)', iostat = status) line

       if ( status .lt. 0 ) exit
       
       if ( line(64:68) .eq. 'gamma' ) then
        write(*,'(A)') "This program is not valid for the gamma version &
                        of VASP."
        write(*,'(A)') "Please rerun with the std version."
        oerror = 2
        return
       endif 

       if ( line(4:11) .eq. 'k-points' ) then
        read(line(30:36),'(I7)') nkpts
        read(line(102:108),'(I7)') nbands
        allocate(kpt(nkpts,3))
        write(*,'(A,I8)') " NBANDS =", nbands
        write(*,'(A,I8)') " NKPTS =", nkpts
       endif 

       if ( line(4:23) .eq. 'dimension x,y,z NGX '.AND. .not.lngx ) then
        lngx = .true.
        read(line(26:30),'(I6)') ngx  
        read(line(37:41),'(I6)') ngy  
        read(line(48:52),'(I6)') ngz  
        allocate(rangex(ngx))
        allocate(rangey(ngy))
        allocate(rangez(ngz))
        write(*,'(A,3I6)') " NGX,Y,Z =", ngx, ngy, ngz
       endif 

       if ( line(4:8) .eq. 'ENCUT' ) then
        read(line(12:18),'(F7.1)') enmax
        write(*,'(A,F7.1)') " ENCUT =", enmax
        if((enmax-enmaxw).gt.0.01) then
         write(*,'(A)') &
          "Cutoff from WAVECAR do not match cutoff from OUTCAR"
          oerror = 7
        elseif(enmax.ne.enmaxw) then
         write(*,'(A)') &
          "Cutoff from WAVECAR do not match cutoff from OUTCAR"
         write(*,'(A,F11.5)') &
          " We will use the value from WAVECAR = ", enmaxw
          enmax = enmaxw
        endif
       endif 

       if ( line(4:8) .eq. 'ISPIN' ) then
        read(line(12:18),'(I7)') nspin
        write(*,'(A,I2)') " NSPIN =", nspin
       endif 

       if ( line(4:7) .eq. 'ISYM' ) then
        read(line(12:18),'(I7)') isym 
        if(isym.ne.0) then
         write(*,'(A,I7)') &
           "ISYM must be equal to 0 and we found ISYM = ", isym
         write(*,'(A)') " Please run VASP again with ISYM = 0."
         oerror = 3
         return
        endif
       endif 

       if ( line(46:77) .eq. 'reciprocal lattice vectors' &
                .AND. .not.cell) then
        do i = 1,3
         read(29,'(45x,3F13.9)') (cell_r(i,j),j=1,3)
        enddo
        cell = .true.
        write(*,'(A)') " Reciprocal unit cell: "
        do i = 1,3
         write(*,'(3F13.9)') (cell_r(i,j),j=1,3)
        enddo 
       endif

       if ( line(2:23) .eq. 'k-points in reciprocal' ) then
        do ik = 1,nkpts
         read(29,'(1x,3F12.8)') (kpt(ik,j),j=1,3)
        enddo
        write(*,'(A)') " K-points in reciprocal lattice:"
        do ik = 1,nkpts
         write(*,'(1x,3F12.8)') kpt(ik,:)
        enddo
       endif

       if ( line(2:8) .eq. 'E-fermi' ) then
        lfermi = .true.
        read(line(11:19),'(F9.4)') efermi
!        write(*,'(A,F9.4)') "E-fermi = ", efermi
       endif 

      enddo
      close(29)

      if(.not.cell) then
       write(*,'(A)') "Could not find the unit cell in OUTCAR."
       oerror = 4
       return
      endif 

      if(.not.lngx) then
       write(*,'(A)') "Could not find NGX, NGY, NGZ in OUTCAR."
       oerror = 5
       return
      endif 

      if(.not.lfermi) then
       write(*,'(A)') "E-fermi not found in OUTCAR. Please check."
       oerror = 6
       return
      endif 

      allocate(onpl(nspin, nkpts))
      onpl = 0

!      enmax = 400.0d0
!      efermi = -4.7847d0
!      nspin = 2
!      nkpts = 2
!      nbands = 47
!      ngx = 98 
!     ngy = 80 
!     ngz = 50

!     cell_r(1,1) = 0.050335570d0
!     cell_r(1,2) = -0.003355705d0
!     cell_r(1,3) = 0.0d0
!     cell_r(2,1) = -0.006711409d0
!     cell_r(2,2) = 0.067114094d0
!     cell_r(2,3) = 0.0d0
!     cell_r(3,1) = 0.0d0
!     cell_r(3,2) = 0.0d0
!     cell_r(3,3) = 0.100000000d0

      
!      kpt = 0
!       kpt(1,:) = [0.01250000,  0.01666667,  0.00000000]
!       kpt(2,:) = [-0.01250000,  0.01666667,  0.00000000]
!       kpt(1,:) = [ 0.25000000,  0.25000000,  0.00000000]
!       kpt(2,:) = [-0.25000000,  0.25000000,  0.00000000]

!      write(*,*) (kpt(1,j),j=1,3)
!      write(*,*) kpt(2,:)

! Generate MappingsCAR      

      call inilpc(ngx, ngy, ngz, rangex, rangey, rangez)

! We write the data to a file to keep compatibility with previous
! versions. The file will be deleted by default.

      open(21,file="MappingsCAR_gen",form="unformatted")

      write(21) efermi
      write(21) nspin, nkpts, ngx, ngy, ngz

      onpl = 0

      spins: do is = 1,nspin
      kpoints: do ik = 1,nkpts
       do iy = 1,ngy
        do ix = 1,ngx
          g1 = rangex(ix) + kpt(ik,1)
          g2 = rangey(iy) + kpt(ik,2)
          
          gix = ( g1*cell_r(1,1) + g2*cell_r(2,1) ) * tpi
          giy = ( g1*cell_r(1,2) + g2*cell_r(2,2) ) * tpi
         
          write(21) gix, giy
        enddo
       enddo 

       do iz = 1, ngz
        do iy = 1, ngy
         do ix = 1, ngx
          g1 = rangex(ix) + kpt(ik,1)
          g2 = rangey(iy) + kpt(ik,2)
          g3 = rangez(iz) + kpt(ik,3)
             
          gix = ( g1*cell_r(1,1) + g2*cell_r(2,1) + g3*cell_r(3,1) ) * tpi
          giy = ( g1*cell_r(1,2) + g2*cell_r(2,2) + g3*cell_r(3,2) ) * tpi
          giz = ( g1*cell_r(1,3) + g2*cell_r(2,3) + g3*cell_r(3,3) ) * tpi
          etotal = (gix**2 + giy**2 + giz**2)/c
          if(etotal.lt.enmax) then
           onpl(is, ik) = onpl(is, ik) + 1
!                  write(17,*) etotal
           write(21) ix,iy,iz
          endif
         enddo
        enddo
       enddo 

      enddo kpoints
      enddo spins

      close(21)
      return
      
      end subroutine mappings_from_outcar

!
! assign ranges
!
      subroutine inilpc (ngx,ngy,ngz,rangex,rangey,rangez)
      implicit none
      integer ngx, ngy, ngz
      integer rangeX(ngx),rangeY(ngy),rangez(ngz)
      integer ix, iy, iz

      do ix = 1,(ngx/2)+1
        rangex(ix) = ix-1
      enddo

      do ix = (ngx/2)+2,ngx
        rangex(ix) = ix-1-ngx
      enddo

      do iy = 1,(ngy/2)+1
        rangey(iy) = iy-1
      enddo

      do iy = (ngy/2)+2,ngy
        rangey(iy) = iy-1-ngy
      enddo

      do iz = 1,(ngz/2)+1
        rangez(iz) = iz-1
      enddo

      do iz = (ngz/2)+2,ngz
        rangez(iz) = iz-1-ngz
      enddo

      return

      end subroutine inilpc

      end module mappingscar_gen
