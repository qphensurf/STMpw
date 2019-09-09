  module currentBardeen
 
CONTAINS

   subroutine matrix_element ( A_S, C_T, ngx, ngy, phi, gx, gy, z, currentSQ)
   use Fourier
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   real (q), parameter :: pi_d = 3.141592653589793238462643383279502884197
   integer, intent (in) :: ngx,ngy
   real (q), intent (in) :: phi,z
   real (q), intent (in) :: gx(ngx*ngy), gy(ngx*ngy)
   complex (q) :: A_S (0:ngx-1,0:ngy-1)
   complex (q) :: C_T (0:ngx-1,0:ngy-1)
   real (q), intent (out) :: currentSQ (0:ngx-1,0:ngy-1)
   complex (q) ::  matrix (0:ngx-1,0:ngy-1)
   real (q) :: decay (0:ngx-1,0:ngy-1)
   integer :: i1, i2, i


! generate all reciprocal lattice vectors:

   i = 0
   do i2 = 0, ngy-1
   do i1 = 0, ngx-1
   i = i+1

! correct decay exponent as a matrix, mapping the
! index with the corresponding reciprocal vector plus k_point

   decay (i1,i2) = sqrt ( 2 * phi + gx(i)**2 + gy(i)**2  )

   enddo
   enddo

   matrix = (0,0)
   matrix(0:ngx-1,0:ngy-1) =  conjg ( C_T(0:ngx-1,0:ngy-1) ) * A_S(0:ngx-1,0:ngy-1) &
           * decay(0:ngx-1,0:ngy-1) * exp(- decay(0:ngx-1,0:ngy-1) * z)


   call fft2d_r (matrix,ngx,ngy)

   currentSQ = real(matrix * conjg (matrix))

   return
   end subroutine matrix_element
!
!
!Tersoff Tip and sample
!
!
   subroutine matrix_TH (A_S,C_T,ngx,ngy,phi,gx,gy,z, TH_t, TH_s)
   use Fourier
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   real (q), parameter :: pi_d = 3.141592653589793238462643383279502884197
   integer, intent (in) :: ngx,ngy
   real (q), intent (in) :: phi,z
   real (q), intent (in) :: gx(ngx*ngy)
   real (q), intent (in) :: gy(ngx*ngy)
   complex (q), intent (in) :: A_S (0:ngx-1,0:ngy-1)
   complex (q), intent (in) :: C_T (0:ngx-1,0:ngy-1)
   real (q), intent (out) :: TH_t (0:ngx-1,0:ngy-1)
   real (q), intent (out) :: TH_s (0:ngx-1,0:ngy-1)
   complex (q) ::  matrix (0:ngx-1,0:ngy-1)
   real (q) :: decay (0:ngx-1,0:ngy-1)
   real (q) :: G (3)
   integer :: i1, i2, i, ix_i, iy_i, ix_f, iy_f


! generate all reciprocal lattice vectors:

   i=0
   do i2 = 0, ngy-1
   do i1 = 0, ngx-1
   i=i+1

! correct decay exponent as a matrix, mapping the
! index with the corresponding reciprocal vector plus k_point

   decay (i1,i2) = sqrt ( 2 * phi + gx(i)**2 + gy(i)**2  )

   enddo
   enddo

!oversized matrix in G
   matrix = (0,0)
   matrix(0:ngx-1,0:ngy-1)= C_T(0:ngx-1,0:ngy-1)   &
          &  * exp(- decay(0:ngx-1,0:ngy-1) * (z))

 
   call fft2d_r (matrix,ngx,ngy)

   TH_t = real(matrix* conjg (matrix))
   matrix = (0,0)

   matrix(0:ngx-1,0:ngy-1) =  A_S * exp(- decay(0:ngx-1,0:ngy-1) * (z))
   matrix = (0,0)
   matrix(0:ngx-1,0:ngy-1)= A_S(0:ngx-1,0:ngy-1)   &
          &  * exp(- decay(0:ngx-1,0:ngy-1) * (z))

   call fft2d_r (matrix,ngx,ngy)

   TH_s = real(matrix * conjg (matrix))


   return
   end subroutine matrix_TH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fermi occupation factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       function Fermi (e)
       implicit none
       integer, parameter :: q = SELECTED_REAL_KIND(10)
       real(q) Fermi, e

        if (e < 0) then
       Fermi= 1.
        else
       Fermi= 0.
        endif

       return
       end function Fermi


  
  end module currentBardeen
