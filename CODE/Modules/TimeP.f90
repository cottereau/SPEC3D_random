module stimeparam 

! Modified by Gaetano 10/02/2005
! Modified by Regis 29/04/2007

type :: time

logical :: acceleration_scheme, velocity_scheme
integer :: ntimeMax, NtimeMin, nSnap, ntrace, ncheck, nenergy,ntime2reverse
real :: alpha, beta, gamma, duration, Time_snapshots_tic,Time_snapshots_toc,Time_snapshots_step, dtmin, rtime, time_energy

end type

end module stimeparam
