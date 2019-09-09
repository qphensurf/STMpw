  module Fourier
!
! In order to make these subroutines work
! please link with lfftw
!
CONTAINS

   subroutine matching (CW, A_G, C_G, ngx, ngy, ngz, z_s, z_t,IGX,IGY,IGZ,NPL,temp)  

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine performs the matching of the coefficients
! A_G at z_s of the FFT(CW)
! C_G at z_t of the FFT(CW)
! A_G and C_G are returned in reciprocal space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Use volumen
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer, parameter :: qs = SELECTED_REAL_KIND(5)
   integer, intent (in):: ngx, ngy, ngz, NPL
   integer :: IGX (NPL), IGY (NPL), IGZ (NPL)
   real (q), intent (in) :: z_s, z_t
   complex (qs), intent (in) :: CW(ngx*ngy*ngz)
   complex (q), intent (out) :: A_G(0:ngx-1,0:ngy-1)
   complex (q), intent (out) :: C_G(0:ngx-1,0:ngy-1)
   complex (q) :: temp(0:ngx-1,0:ngy-1,0:ngz-1),caca
   integer :: ix, iy, iz ,l, i_s, i_t,test

   temp = (0,0)
   A_G = (0,0)
   C_G = (0,0)
   do l = 1, NPL
   temp (IGX(l)-1,IGY(l)-1,IGZ(l)-1) = CW(l)*1.0_q
   enddo

   call fft1d_r (temp, ngx, ngy, ngz)



! find index for the matching

   i_s = mod (int(z_s*ngz), ngz)
   i_t = mod (int(z_t*ngz), ngz)


! matching at the corresponding plane

   A_G(0:ngx-1,0:ngy-1) = temp (0:ngx-1,0:ngy-1,i_s-1) ! sample
   C_G(0:ngx-1,0:ngy-1) = temp (0:ngx-1,0:ngy-1,i_t-1) ! tip


   return
   end subroutine matching 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft1d_r (temp, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer i_x, i_y,ngx, ngy, ngz
   complex (q) matrix (0:ngz-1)
   complex (q) temp (0:ngx-1,0:ngy-1,0:ngz-1)

   do i_x = 0, ngx-1
   do i_y = 0, ngy-1
   matrix(:) = temp (i_x,i_y,:)
   call dfftw_plan_dft_1d (plan,ngz,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   temp (i_x,i_y,:) = matrix(:)

   enddo
   enddo

   end subroutine fft1d_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3-D fft FORWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft3d (matrix, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy, ngz
   complex (q) matrix (0:ngx-1,0:ngy-1,0:ngz-1)

   call dfftw_plan_dft_3d (plan,ngx,ngy,ngz,matrix,matrix,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft3d 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft3d_r (matrix, ngx, ngy, ngz)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy, ngz
   complex (q) matrix (0:ngx-1,0:ngy-1,0:ngz-1)

   call dfftw_plan_dft_3d (plan,ngx,ngy,ngz,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft3d_r
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2-D fft FORWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft2d (matrix, ngx, ngy)
   use declfft
   implicit none
   integer ngx, ngy
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   complex (q) matrix (0:ngx-1,0:ngy-1)

   call dfftw_plan_dft_2d (plan,ngx,ngy,matrix,matrix,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)

   end subroutine fft2d 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2-D fft BACKWARD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft2d_r (matrix, ngx, ngy)
   use declfft
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer ngx, ngy
   complex (q) matrix (0:ngx-1,0:ngy-1)

   call dfftw_plan_dft_2d (plan,ngx,ngy,matrix,matrix,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan)
        call dfftw_destroy_plan(plan)
   end subroutine fft2d_r 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine density (density1,CW,ngx, ngy, ngz, IGX,IGY,IGZ,NPL)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer, parameter :: qs = SELECTED_REAL_KIND(5)
   integer, intent (in):: ngx, ngy, ngz, NPL
   integer :: IGX (NPL), IGY (NPL), IGZ (NPL)
   complex (qs), intent (in) :: CW(ngx*ngy*ngz)
   real (q), intent (out) :: density1(ngz)
   complex (q) :: temp(0:ngx-1,0:ngy-1,0:ngz-1)
   integer :: ix, iy, iz ,l, i_s, i_t,test

   temp = (0,0)
   density1 = 0
   do l = 1, NPL
   temp (IGX(l)-1,IGY(l)-1,IGZ(l)-1) = CW(l)*1.0_q
   enddo


   call fft3d_r (temp, ngx, ngy, ngz)

   do iz =1, ngz
   do ix = 1, ngx
   do iy = 1, ngy
   density1(iz) = density1(iz) + real(temp (ix-1,iy-1,iz-1)*conjg(temp (ix-1,iy-1,iz-1)))
   enddo
   enddo
   enddo

   density1 = density1/ (ngx*ngy)
   end subroutine density

  end module Fourier
