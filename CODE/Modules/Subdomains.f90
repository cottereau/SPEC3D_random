 module ssubdomains
! Gaetano Festa 01/02/2005

type Subdomain

logical :: Filtering, Px, Py, Pz, Left, Forward, Down

integer :: NGLLx, NGLLy, NGLLz, wpml, npow

real :: Pspeed, Sspeed, Ddensity, Dt, Apow, freq, DLambda, DMu,coPMLrate

real,dimension(21)::Coefs21

real, dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
real, dimension (:,:), pointer :: hprimex, hTprimex
real, dimension (:), pointer :: GLLcy, GLLpoly, GLLwy
real, dimension (:,:), pointer :: hprimey, hTprimey
real, dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
real, dimension (:,:), pointer :: hprimez, hTprimez

character(len=2) :: material_type

logical :: random,aniso, reload_mat
character(len=1) :: random_type, Correlation_type,random_type_kernel
real :: CorLx, CorLy, CorLz, delta,delta_kernel
integer :: ChosenSeed, CorNx, CorNy, CorNz
integer :: PML_GLL_out=0

logical :: SpheDev
!real,dimension(:,:),allocatable::Spheric,Deviatoric


end type

contains

subroutine Lame_coefficients (S)

type (Subdomain) :: S

S%DMu = S%Sspeed**2 * S%Ddensity
S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed**2 ) * S%Ddensity 

!!$if (S%SpheDev) then
!!$   allocate(S%Spheric(1:6,1:6))
!!$   allocate(S%Deviatoric(1:6,1:6))
!!$   S%Spheric=0.
!!$   S%Spheric(1:3,1:3)=1.
!!$   S%Spheric(:,:)=1./3.*S%Spheric
!!$   S%Deviatoric(1:6,1:6)=1.
!!$   S%Deviatoric=S%Deviatoric-S%Spheric
!!$   S%Deviatoric(4:6,4:6)=S%Deviatoric(4:6,4:6)/sqrt(2.)
!!$endif

end subroutine

end module ssubdomains
