  module declaration

  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: qs = SELECTED_REAL_KIND(5)

! sqrt(pi):
  real (q), parameter :: sqrpi = 1.7724538509055151137
  real (q), parameter :: pi_d = 3.141592653589793238462643383279502884197
  real (q), parameter :: bohr = 0.529177210903
  real (q), parameter :: hartree = 27.211386245988


!SCALARS:

! integers 
  integer :: Ngx, Ngy, Ngz, n_species
  integer :: Number_of_SPIN, Total_Bands, Number_of_KPOINTS
  integer ::  unitWF, unitI, unitdI,unitMAP,unitTH,unitTIP,unitTHcub,unitdITHcub,unitIV
  integer ::  unitIdat, unitdIdat, unitTHdat, unitdITHdat, unitTIPdat
  integer ::  spin1,  kp, iband, tband, sband
  integer ::  ios, N_sampling_z 
  integer :: ix, iy, iz, ip
  integer :: NSPIN, NKPTS, spin_factor, nat, nV, jv, ndiv, npts, oerror
! From VASP
      INTEGER, PARAMETER :: ICMPLX=16,MRECL=10000000
      INTEGER IDUM, I, J, ISP, ISPIN
      INTEGER NPL, IREC
      INTEGER ISPINOR, ISPINOR_, NRSPINORS, NPRO1, NPRO1_
      INTEGER NIS, NT, NTYP, LMMAXC, NI, NPRO
      INTEGER LMDIMp,NIONSp,NCDIJp
      INTEGER LMDIMf,NIONSf,NCDIJf
      INTEGER LMDIM,NIONS,NCDIJ
      INTEGER NPRODp, NPROp, NTYPp,NPROD
      INTEGER NPRODf, NPROf, NTYPf

! reals 
  real (q) :: Efermi,phi,Zmax,Zsurf,Zmin,z
  real (q) :: z_s, z_t, Fermi_t, Fermi_s
  real (q) :: sigma = 0.07/27.2116
  real (q) :: stepX, stepY, stepZ, Ztip,offset,distance
  real (q) :: constant_I, constant_dIdV
  real (q) :: delta_E,vol,surf,CutoffEnergy
  real (q) :: lat_par, cell(3,3)
  real (q) :: start, finish, t_wf, st_tmp, t_ma, fi_tmp
  real (q) :: Vmin, Vmax, factor


! From VASP
      REAL(q)  RDUM, RISPIN, ENMAXF
      REAL(q)  RTAG, RKPTSF, RBANDF, RNPL
      REAL(q)  RKPTS, RBANDS, ENMAXW
      REAL(q) A(3,3), UCELL(3,3)
      real(q) vec(3,3)

!logical
  LOGICAL :: LGAMMA, fichero, Bardeen, LDIDV, MAPfile
  LOGICAL :: wsxm = .false., dat = .false., cube = .false.

! characters
 character ( len = 99) :: NameWF,NameMAP,NamePOSCAR,name_file,cmd
 character ( len = 6) :: volt,cnpt



!ARRAYS:

! integers
  integer, allocatable :: iGX(:), iGY(:), iGZ(:), zat(:), ngp(:,:) 
  integer, pointer :: onpl(:,:)
! From VASP
      INTEGER, ALLOCATABLE:: LMMAX (:), NITYP (:)
      INTEGER, ALLOCATABLE:: LMMAXp (:), NITYPp (:)
      INTEGER, ALLOCATABLE:: LMMAXf (:), NITYPf (:)


! reals 
  real (q), allocatable :: CQIJ (:,:,:,:)
  real (q), allocatable :: gx(:), gy(:), gz(:), cngx(:,:), cngy(:,:)
  real (q), allocatable :: density_z(:), density1(:)
  real (q), allocatable :: coord(:,:), V(:), V1(:), IV(:,:), cdIdV(:,:)

! From VASP
      REAL(q), ALLOCATABLE:: VKPT(:,:)
      REAL(q), ALLOCATABLE :: FERTOT(:,:,:)
      REAL(q), ALLOCATABLE :: currentSQ (:,:), intensity (:,:,:,:)
      REAL(q), ALLOCATABLE :: intensity2 (:,:)
      REAL(q), ALLOCATABLE :: Tersoff_s (:,:,:,:)
      REAL(q), ALLOCATABLE :: Tersoff_s2 (:,:)
      REAL(q), ALLOCATABLE :: Tersoff_t (:,:,:,:)
      REAL(q), ALLOCATABLE :: Tersoff_c (:,:,:,:)
      REAL(q), ALLOCATABLE :: TH_t (:,:)
      REAL(q), ALLOCATABLE :: TH_s (:,:)
      REAL(q), ALLOCATABLE :: dintensity_dV (:,:,:,:)


! complex 
  complex (q), allocatable :: EIG(:)
  complex (qs), allocatable :: CW(:)
  complex (q), allocatable :: A_G (:,:,:), C_G (:,:,:)
  complex (q), allocatable ::  A_S (:,:), C_T (:,:), temp(:,:,:)


  end module declaration
