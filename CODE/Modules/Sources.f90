module ssources

! Modified by Gaetano Festa 31/01/2005
! Modified by Paul Cupillard 05/06/2005


type :: elem_source
!!$integer :: nr !number of elements containning the source either in its interior or on its boundary
double precision :: eta,xi,zeta
!!$real, dimension (0:2,0:2) :: Scoeff
real, dimension (:,:,:), pointer :: ExtForce
real, dimension (:,:,:,:), pointer :: explosion
integer,pointer::indspecel
end type 

type :: proc_source
logical,pointer::source_affected
integer,pointer :: nr !number of elements containning the source either in its interior or on its boundary
type(elem_source),dimension(:),pointer::elem
endtype

type :: Source
integer :: i_dir, i_type_source, i_time_function
real :: Xsource,Zsource,Ysource, tau_b,cutoff_freq
!!$logical,dimension(:),pointer::proc
!!$type(elem_source), dimension(:), pointer :: Elem
type(proc_source),dimension(:),pointer::proc
end type 

!!$type :: Source
!!$
!!$integer :: i_dir, i_type_source, i_time_function, proc, elem
!!$integer , dimension(0:2) :: gll
!!$real :: Xsource,YSource,Zsource, tau_b,cutoff_freq
!!$
!!$
!!$end type 

contains

! ###############################################
  real function CompSource (Sour,time,np)

    type (source) :: Sour
    integer ::np
    real :: time

!!$print *, Sour%i_type_source, 'Sour%i_type_source'
    select case (Sour%i_type_source)
    case (1)   ! Impulse
       if (np /= Sour%i_dir) then
          CompSource = 0. 
       else   
          select case (Sour%i_time_function)
          case (1) 
             CompSource = Gaussian (time,Sour%tau_b)
          case (2)
             CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
          end select
       end if
    case (2) ! Explosion
       select case (Sour%i_time_function)
       case (1)
          CompSource = Gaussian (time,Sour%tau_b)
       case (2)
          CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
       end select
    end select
    return
  end function CompSource

! ################################################
real function Gaussian (time, tau)

real :: tau,time

!!$Gaussian = -(time-tau) * exp (-(time-tau)**2/tau**2)

Gaussian = exp (-(time-tau)**2/0.01)
return
end function
 
! ################################################
real function Ricker (time,tau,f0)

use pig

real :: time, tau, f0
real :: sigma

sigma = pi * f0 * (time-tau)
sigma = sigma**2
Ricker = (1 - 2*sigma) * exp(-sigma)

return
end function

! #################################################
end module ssources
