program Imagen_gwy
! Program for preparing input data for Gwyddion
!
! Load Topography.gwy as xyz with 1e-10 as Lateral units 
!  and 1e-10 for Value units for constant current and 1 for constant height.
! Then apply XYZData -> Rasterize with the appropiate Resolution 
!
       implicit none
            
       integer, parameter :: q = SELECTED_REAL_KIND(10)
       integer, parameter :: qs = SELECTED_REAL_KIND(5)
       
       integer :: N_sampling_z, ngy, ngx, iz, iy, ix, i, j
       integer :: var, ix2, idx, idy
       integer :: nx, ny, jx, jy
       integer, parameter :: unit_topo = 1
       real (q) :: plotting_current, Hei
       real (q), allocatable :: intensity (:,:,:)
       real (q), allocatable :: z (:), x (:), y (:)
       real (q), allocatable :: zvalue (:,:)
       real (q), allocatable :: temp (:,:)
       logical :: Bardeen, present, presentgz
       character (len = 80) :: line, name_file
       character (len = 6) :: name
       double precision :: alat, lat(3,3), mdl, xplot, yplot, ymax

! Read unit cell from POSCAR

       line='POSCAR'
       inquire(file=line, exist=present)
       if(present) then
         open(3,file='POSCAR',status='old')
       elseif(.NOT.present) then
         line='../POSCAR'
         inquire(file=line, exist=present)
         if(present) then
           open(3,file='../POSCAR',status='old')
         elseif(.NOT.present) then
           write(*,'(A21)') 'POSCAR is not present'
           stop
         endif  
       endif

       read(3,'(A80)') line
       read(3,'(F19.14)') alat
       do i=1,3
        read(3,*) (lat(i,j),j=1,3)
       enddo

       mdl = dsqrt(lat(1,1)**2+lat(1,2)**2)
       lat(1,1) = lat(1,1)/mdl 
       lat(1,2) = lat(1,2)/mdl 

       mdl = dsqrt(lat(2,1)**2+lat(2,2)**2)
       lat(2,1) = lat(2,1)/mdl 
       lat(2,2) = lat(2,2)/mdl 
       close(3)
!
       print *, 'T for Bardeen and F for Tersoff-Hamann: '
       read (*,*) Bardeen

       if(.not.Bardeen) then
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
       endif  

       if(Bardeen) then
         line='Bardeen.dat'
         inquire(file=line, exist=present)
         if(.NOT.present) then
           write(*,'(A11,A15)') line,' is not present'
           stop
         endif
         open(2,file='Bardeen.dat',status='old')
       endif  

       read (2,*) N_sampling_z, ngy, ngx

       allocate (intensity (0:ngx-1, 0:ngy-1, N_sampling_z))
       allocate (z(N_sampling_z))
       allocate (y(0:ngy))
       allocate (x(0:ngx))
       allocate (zvalue(0:ngx, 0:ngy))
       allocate (temp(0:2*ngx, 0:2*ngy))

       do iz= 1, N_sampling_z
         do iy= 0, ngy-1
           do ix = 0, ngx-1
             read (2, '(4g14.4)') z(iz), y(iy), x(ix), intensity(ix,iy,iz)
           enddo
         enddo
       enddo

       close (2)
       if(presentgz) call system ("gzip TH.dat")

! test of ordering of coordinate
       if (z(1) > z(2)) then
         print *, 'Please change test in currents'
         print *, 'in this code, because we have'
         print *, 'assumed that z increases in STM.dato'
         stop
       endif

       print *, 'Number of repetitions in x ?='
       read (*,*) nx

       print *, 'Number of repetitions in y ?='
       read (*,*) ny
         
       print *, 'either constant height (type 1) OR'
       print *, 'either constant current (type 2) '
       read (*,*) var

       if(var.eq.1) then
        print *, 'Value of the plotting HEIGHT ?='
        read (*,*) Hei

        do ix = 0, ngx-1
          do iy= 0, ngy-1
            do iz= 1, N_sampling_z-1
               if ( Hei >= z(iz) .and. &
         &          Hei <= z(iz+1) ) then

           zvalue (ix, iy) =  intensity(ix,iy,iz) + &
         &  (Hei-z(iz)) * (intensity(ix,iy,iz+1)-intensity(ix,iy,iz))/&
         &  (z(iz+1)-z(iz))

                endif
            enddo
          enddo
        enddo
       elseif(var.eq.2) then

        print *, 'Value of the plotting current ?='
        read (*,*) plotting_current
        write (*,*) plotting_current

        do ix = 0, ngx-1
          do iy= 0, ngy-1
            do iz= 1, N_sampling_z-1
              if ( plotting_current <= intensity(ix,iy,iz) .and. &
                  plotting_current >= intensity(ix,iy,iz+1) ) then
                zvalue (ix, iy ) = z (iz) + &
                 (plotting_current-intensity(ix,iy,iz))*(z(iz+1) - z(iz))/ &
                 (intensity(ix,iy,iz+1)-intensity(ix,iy,iz))
              endif
            enddo
          enddo
        enddo
       endif

! printing out results

       y(ngy) = y(ngy-1) + (ny-1)*(y(ngy-1)+y(1))+y(1)
       ymax = (x(ngx-1) + (nx-1)*(x(ngx-1)+x(1))+x(1))*lat(1,2) + y(ngy)*lat(2,2)

       if(Bardeen) then
         if(var.eq.1) then
           write(name,'(F4.1)') Hei
           name_file = 'Topography_Bardeen_ch_'//trim(adjustl(name))//'A.gwy'      
         else
           write(name,'(ES6.0)') plotting_current
           name_file = 'Topography_Bardeen_cc_'//trim(adjustl(name))//'.gwy'      
         endif
       else
         if(var.eq.1) then
           write(name,'(F4.1)') Hei
           name_file = 'Topography_ch_'//trim(adjustl(name))//'A.gwy'      
         else
           write(name,'(ES6.0)') plotting_current
           name_file = 'Topography_cc_'//trim(adjustl(name))//'.gwy'      
         endif
       endif

       open (unit_topo, file = name_file)

       do ix = 1, nx
        do jx= 0, ngx-1
         do iy = 1, ny
          do jy= 0, ngy-1
           xplot = x(jx) + (ix-1)*(x(ngx-1)+x(1))
           yplot = y(ngy) - y(jy) - (iy-1)*(y(ngy-1)+y(1))
           yplot = y(jy) + (iy-1)*(y(ngy-1)+y(1))
           write (unit_topo, *) xplot*lat(1,1) + yplot*lat(2,1), &
                         ymax - xplot*lat(1,2) - yplot*lat(2,2), zvalue (jx,jy)
          enddo
         enddo
         xplot = x(jx) + (ix-1)*(x(ngx-1)+x(1))
         yplot = y(ngy)
         write (unit_topo, *) xplot*lat(1,1) + yplot*lat(2,1), &
                         ymax - xplot*lat(1,2) - yplot*lat(2,2), zvalue (jx,0)
         write (unit_topo, *) 
        enddo
       enddo
       do iy = 1, ny
        do jy= 0, ngy-1
         xplot = x(ngx-1) + (nx-1)*(x(ngx-1)+x(1))+x(1)
         yplot = y(ngy) - y(jy) - (iy-1)*(y(ngy-1)+y(1))
         yplot = y(jy) + (iy-1)*(y(ngy-1)+y(1))
         write (unit_topo, *) xplot*lat(1,1) + yplot*lat(2,1), &
                         ymax - xplot*lat(1,2) - yplot*lat(2,2), zvalue (0,jy)
        enddo
       enddo
       xplot = x(ngx-1) + (nx-1)*(x(ngx-1)+x(1))+x(1)
       yplot = y(ngy)
       write (unit_topo, *) xplot*lat(1,1) + yplot*lat(2,1), &
                         ymax - xplot*lat(1,2) - yplot*lat(2,2), zvalue (0,0)
       close (unit_topo)

     stop
end program Imagen_gwy
