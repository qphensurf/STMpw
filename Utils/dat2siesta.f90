program dat2siesta
! Program to transform .dat in .siesta files
!
       implicit none
            
       integer, parameter :: q = SELECTED_REAL_KIND(10)
       
       real (q), parameter :: bohr = 0.529177210903

       integer :: N_sampling_z, ngy, ngx, iz, iy, ix, i, j
       real (q), allocatable :: intensity (:,:,:)
       real (q), allocatable :: z (:), x (:), y (:)
       real (q) :: factor, Zmax
       real (q) :: alat, lat(3,3)
       logical :: present, presentgz
       character (len = 80) :: line

       ! Read unit cell from POSCAR

       line='POSCAR'
       inquire(file=line, exist=present)
       if(.NOT.present) then
         write(*,'(A6,A15)') line,' is not present'
       stop
       endif

       open(3,file='POSCAR',status='old')
       read(3,'(A80)') line
       read(3,'(F19.14)') alat
       do i=1,3
        read(3,*) (lat(i,j),j=1,3)
       enddo
       lat = lat * alat / bohr

       close(3)

       line='TH.dat'
       inquire(file=line, exist=present)
       if(present) then
         open(2,file='TH.dat',status='old')
       elseif(.NOT.present) then  
         line='TH.dat.gz'
         inquire(file=line, exist=presentgz)
         if(presentgz) then
           call system ("gunzip TH.dat.gz")
           open(2,file='TH.dat',status='old')
         else  
           write(*,'(A6,A15)') line,' is not present'
           stop
         endif
       endif  
       
       write(*,*) "Multiplying factor for WSxM ?"
       read(*,*) factor

       write(*,*) "Zmax, maximum tip-surface distance (from 'Zsurf') in angstroms?"
       read(*,*) Zmax  
       Zmax = Zmax / bohr
       lat(3,3) = Zmax

       open(3,file='TH.siesta',form = 'unformatted')
       
       write (3) lat  

       read (2,*) N_sampling_z, ngy, ngx
       
       write(3) ngx, ngy, N_sampling_z,1

       allocate (intensity (0:ngx-1, 0:ngy-1, N_sampling_z))
       allocate (z(N_sampling_z))
       allocate (y(0:ngy))
       allocate (x(0:ngx))

       do iz= 1, N_sampling_z
         do iy= 0, ngy-1
           do ix = 0, ngx-1
             read (2, '(4g14.4)') z(iz), y(iy), x(ix), intensity(ix,iy,iz)
           enddo
           write (3) (real(intensity(ix,iy,iz)*factor), ix=0,ngx-1)
         enddo
       enddo

       close (2)
       if(presentgz) call system ("gzip TH.dat")

     stop
end program dat2siesta
