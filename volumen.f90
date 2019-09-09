! module that contains subroutines for
! performing vectorial products and
! volume evaluations
  module volumen

CONTAINS
   subroutine surface (vec1,vec2,surf)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   real (q) :: vec1(3)
   real (q) :: vec2(3)
   real (q) ::  product(3)
   real (q) :: surf

   call vectorial_product (vec1,vec2,product)

   surf = sqrt( sum ( product (:)**2 ) )

   return
   end subroutine surface


   subroutine volume (vec1,vec2,vec3,vol)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   real (q) :: vec1(3)
   real (q) :: vec2(3)
   real (q) :: vec3(3), product(3)
   real (q) :: vol

   call vectorial_product (vec1,vec2,product)

   vol = ( sum ( vec3 (:) * product (:) ) )

   return
   end subroutine volume

   subroutine vectorial_product (vec1,vec2,product)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   real (q) :: vec1(3)
   real (q) :: vec2(3)
   real (q) :: product(3)

   product (1) = vec1(2) * vec2(3) -  &
             & vec1(3) * vec2(2)
   product (2) = vec1(3) * vec2(1) -  &
             & vec1(1) * vec2(3)
   product (3) = vec1(1) * vec2(2) -  &
             & vec1(2) * vec2(1)

   return
   end subroutine vectorial_product 

        subroutine  reciprocal(a1,a2,a3,vol,vec1,vec2,vec3)
        implicit none
        integer, parameter :: q = SELECTED_REAL_KIND(10)
        real (q), parameter :: PI_D = 3.141592653589793238462643383279502884197
        real (q) vol
        real (q) pv(3),a1(3),a2(3),a3(3)
        real (q) vec1(3),vec2(3),vec3(3)

        call vectorial_product(a2,a3,pv)
              vec1=2*PI_D*pv/vol
        call vectorial_product(a3,a1,pv)
              vec2=2*PI_D*pv/vol
        call vectorial_product(a1,a2,pv)
              vec3=2*PI_D*pv/vol

        return
        end subroutine  reciprocal
        

   subroutine Cutoff_sphere (ix,iy,iz,CutoffEnergy,VKPT,vec,test,ngx,ngy,ngz)
   implicit none
   integer, parameter :: q = SELECTED_REAL_KIND(10)
   integer :: ix,iy,iz,test
   real (q), intent (in):: CutoffEnergy
   real (q) :: Energie_cinetique
   real (q) :: VKPT (3), vec (3,3) !k-points and reciprocal unit cell
   real (q) :: G(3) ! G-vector
   integer :: ngx, ngy, ngz, iGx, iGy, iGz


! initialisation

    test = 0

! vecteur G correspondant
    iGx = mod (ix-1,ngx/2)-(ix-1)/(ngx/2)*(ngx/2)
    iGy = mod (iy-1,ngy/2)-(iy-1)/(ngy/2)*(ngy/2)
    iGz = mod (iz-1,ngz/2)-(iz-1)/(ngz/2)*(ngz/2)

    G(:)=(iGx+VKPT(1))*vec(:,1)+(iGy+VKPT(2))*vec(:,2)+ (iGz+VKPT(3))*vec(:,3)

    Energie_cinetique = 0.5 * sum( G (:)**2 )

    if (Energie_cinetique < CutoffEnergy) then

        test = 1

    endif

   return
   end subroutine  Cutoff_sphere

  end module volumen
