program Imagen_gnu
! program for preparing input data for gnuplot
!
! gnuplot>set pm3d
! gnuplot>splot "Topography.gnu" u 1:2:3 with lines
!
       implicit none
            
       integer, parameter :: q = SELECTED_REAL_KIND(10)
       integer, parameter :: qs = SELECTED_REAL_KIND(5)
       
       integer :: N_sampling_z, ngy, ngx, iz, iy, ix
       integer :: var, ix2, idx, idy
       integer :: nx, ny, jx, jy
       integer, parameter :: unit_topo = 1
       real (q) :: plotting_current, Hei
       real (q), allocatable :: intensity (:,:,:)
       real (q), allocatable :: z (:), x (:), y (:)
       real (q), allocatable :: zvalue (:,:)
       real (q), allocatable :: temp (:,:)
       logical :: Bardeen

       print *, 'T for Bardeen and F for Tersoff-Hamann: '
       read (*,*) Bardeen
       if(.not.Bardeen) open(2,file='TH.dat',status='old')
       if(Bardeen) open(2,file='Bardeen.dat',status='old')

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

!       do ix = 0, ngx-1
!         do iy= 0, ngy-1
!           if (zvalue (ix, iy ) < 5.5 ) then
!             zvalue (ix, iy ) = 5.5
!             print *, 'WARNING: we saturate minimum distance at 5.5 Ang'
!           endif
!         enddo
!       enddo

! Shift to middle of the cell
       var = 0
       print *, ' Type (1) if you want ot shif the cell='
       read (*,*) var  

       if (var == 1) then
         temp = 0
         print *, ' Type (1) if you want to shift the molecule in a1'
         print *, ' Type (2) if you want to shift the molecule in a2'
         print *, ' Type (3) if you want to shift to the middle of the cell'
         read (*,*) var

         select case (var)
           case (1)
                idx = ngx / 2
                idy = 0
           case (2)
                idx = 0
                idy = ngy/2
           case (3)
                idx = ngx / 2
                idy = ngy / 2
         end select
       
         do ix = 0, ngx-1
           do iy = 0, ngy-1
             print *, ix, mod (ix+idx, ngx)
             temp (ix,iy) = zvalue (mod (ix+idx, ngx), mod (iy+idy,ngy))
           enddo
         enddo

         do ix2 = 0, ngx - 1
           do iy = 0, ngy - 1
             zvalue (ix2, iy) = temp (ix2, iy)
           enddo
         enddo

       endif

! printing out results

       open (unit_topo, file = 'Topography.gnu')
       x(ngx) = x(ngx-1)+x(1)
       y(ngy) = y(ngy-1)+y(1)
       zvalue(:,ngy)=zvalue(:,0)
       zvalue(ngx,:)=zvalue(0,:)
       zvalue(ngx,ngy)=zvalue(0,0)
       do ix = 1, nx
        do jx= 0, ngx
         do iy = 1, ny
          do jy= 0, ngy
           write (unit_topo, *) x(jx) + (ix-1)*(x(ngx-1)+x(1)), &
                                y(jy) + (iy-1)*(y(ngy-1)+y(1)), zvalue (jx,jy)
          enddo
         enddo
         write (unit_topo, *) 
        enddo
       enddo
       close (unit_topo)

     stop
end program Imagen_gnu
