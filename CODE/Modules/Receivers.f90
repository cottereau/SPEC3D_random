module sreceivers

type :: receiver

! Modified by Gaetano Festa 31/01/2005 
! Modified by Paul Cupillard 05/06/2005

integer :: elem, proc, flag
integer , dimension(0:2) :: gll
real :: xRec,yRec,zRec, xi,eta,zeta
real, dimension(:,:), allocatable:: StoreTrace
real, dimension(:,:,:), allocatable:: pol
real, dimension(:,:,:,:), allocatable:: coeff

end type 
end module sreceivers
