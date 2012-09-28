program  drive_sem

! modified Regis 04/2007

use sdomain
!use MPI

implicit none

include 'mpif.h'

type(domain), target :: Tdomain

integer :: code, rg, nb_procs, ntime, i_snap, n,i,j, icount, icountc, icountE, N_step_0
character (len=100) :: fnamef
real::T_comp0, T_comp1,T_comp2

call mpi_init(code) 
call mpi_comm_rank(mpi_comm_world, rg, code)
call mpi_comm_size(mpi_comm_world, nb_procs, code)

! ##############  Begin of the program  #######################
write (*,*) "Define mesh properties ",rg
call read_input (Tdomain,rg)

print*,'SO',Tdomain%logicD%super_object_local_present
print*,'Neu',Tdomain%logicD%neumann_local_present

if (Tdomain%logicD%super_object_local_present) then
 if (Tdomain%super_object_type == "P") then
!   write(*,*) "Define Plane Wave properties",rg
   call define_planew_properties (Tdomain)
 endif
endif
if (Tdomain%logicD%neumann_local_present) then
!  write(*,*) "Define Neumann properties",rg
  call define_neu_properties (Tdomain)
endif

!write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes ",rg
call compute_GLL (Tdomain)

!write (*,*) "Define a global numbering for the collocation points ",rg
call global_numbering (Tdomain,rg)

!write (*,*) "Computing shape functions within their derivatives ",rg
if (Tdomain%n_nodes == 8) then
      call shape8(TDomain,rg)   ! Linear interpolation
else if (Tdomain%n_nodes == 27) then
      write (*,*) "Sorry not yet implemented in the code... Wait for an upgrade ",rg
      stop
!      call shape27(TDomain)  ! Quadratic interpolation
else
      write (*,*) "Bad number of nodes for hexaedral shape ",rg
      stop
endif

!write (*,*) "Compute Courant parameter ",rg
call compute_Courant (Tdomain,rg)

if (Tdomain%any_PML)   then
!    write (*,*) "Attribute PML properties ",rg
    call PML_definition (Tdomain,rg) 
endif

if (Tdomain%logicD%any_source) then
   do n=0,Tdomain%n_source-1
      allocate(Tdomain%sSource(n)%proc(0:nb_procs-1))
      do i=0,nb_procs-1
         allocate(Tdomain%sSource(n)%proc(i)%source_affected)
         Tdomain%sSource(n)%proc(i)%source_affected=.false.
      enddo
   enddo
   call mpi_barrier(MPI_COMM_WORLD,code)
   call SourcePosition(Tdomain,rg,nb_procs)
endif

!!$if (Tdomain%logicD%any_source) then
!!$!    write (*,*) , "Computing point-source parameters ",rg
!!$    call SourcePosition(Tdomain,nb_procs)
!!$endif



if (Tdomain%logicD%save_trace) then
!    write (*,*) "Computing receiver locations ",rg
    call ReceiverPosition(Tdomain,rg)
endif

!write (*,*) "Allocate fields ",rg
call allocate_domain (Tdomain,rg)
T_comp0=MPI_WTIME()
!print *, 'DriveSem::Define_array_begin',T_comp,rg
	write (*,*) "Compute mass matrix and internal forces coefficients ",rg
call define_arrays (Tdomain,rg)
	T_comp1=MPI_WTIME()
if (rg==0) print *, 'DriveSem::Generation time',T_comp1-T_comp0,'seconds'

	write (*,*) "FPML",Tdomain%any_FPML,rg

!if ( Tdomain%logicD%run_exec ) then
write (*,*) "Entering the time evolution ",rg
Tdomain%TimeD%rtime = 0
Tdomain%TimeD%NtimeMin = 0

if ((Tdomain%logicD%save_snapshots) .and. (Tdomain%Field_Order(4))) Tdomain%TimeD%time_snapshots_step=Tdomain%TimeD%dtmin

if (Tdomain%logicD%run_restart) call read_restart(Tdomain,rg)

if (Tdomain%logicD%save_snapshots) then
 Tdomain%timeD%nsnap = Tdomain%TimeD%time_snapshots_step / Tdomain%TimeD%dtmin
 N_step_0=Tdomain%TimeD%time_snapshots_tic / Tdomain%TimeD%dtmin 
 icount = 0
endif
if (Tdomain%logicD%save_energy) then
    Tdomain%timeD%nenergy = Tdomain%TimeD%time_energy / Tdomain%TimeD%dtmin
    icountE = 0
endif

icountc = 0
i_snap=1
do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1
   call Newmark (Tdomain, rg, ntime)

   Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

   if ((Tdomain%logicD%save_snapshots) &
        .and. (Tdomain%TimeD%rtime .gt. Tdomain%TimeD%time_snapshots_tic) &
        .and. (Tdomain%TimeD%rtime .lt. Tdomain%TimeD%time_snapshots_toc)) then

      i_snap = mod (ntime - N_step_0, Tdomain%TimeD%nsnap)
!print *, 'drive_Sem::nsnap',Tdomain%TimeD%nsnap
      if (i_snap == 0) then
         icount = icount+1
         call savefield (Tdomain,ntime,rg,icount)
         !            write (*,*) "A new snapshot is recorded ",rg
      endif
   endif

   if (Tdomain%logicD%save_energy) then
      if ((mod (ntime+1, Tdomain%TimeD%nenergy)) == 0) then
         icountE = icountE+1
         call saveenergy (Tdomain,rg,icountE)
         !            write (*,*) "A new energy snapshot is recorded ",rg
      endif
   endif

   if (Tdomain%logicD%save_trace) then
      if (mod (ntime+1,Tdomain%TimeD%ntrace) == 0) then
        call savetrace (Tdomain,rg,int(ntime/Tdomain%TimeD%ntrace))
        do n = 0, Tdomain%n_receivers-1
          if (rg == Tdomain%sReceiver(n)%proc) deallocate (Tdomain%sReceiver(n)%StoreTrace)
        enddo
      endif
    endif

    ! Checkpoint restart
    if (Tdomain%logicD%save_restart)  then
      if ( mod (ntime+1,Tdomain%TimeD%ncheck) == 0 ) then
         icountc = icountc + 1
         call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,ntime+1,rg,icountc)
      endif
    endif
enddo
!print *,'DRIVE_SEM::check2'
!print *,'DRIVE_SEM::Field_Order',Tdomain%Field_Order(4:4),rg
if ((Tdomain%logicD%save_snapshots) .and. (Tdomain%Field_Order(4))) then
   Tdomain%Field_Order(1)=.false.
   Tdomain%Field_Order(2)=.false.
   Tdomain%Field_Order(3)=.false.
   Tdomain%Field_Order(5)=.false.

!print *,'DRIVE_SEM::check2a'
   Tdomain%replay=.true.
!print *,'DRIVE_SEM::check2b'
   !print *,'DRIVE_SEM::check3'
Tdomain%TimeD%rtime=0.
   do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1
      call Newmark (Tdomain, rg, ntime)
	!print *, 'DRIVE_SEM::check4'
      Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

      !if (Tdomain%logicD%save_snapshots) then
         i_snap = mod (ntime+1, Tdomain%TimeD%nsnap) 
         if (i_snap == 0) then
!print *,'drivesem', ntime, Tdomain%TimeD%nsnap
            icount = icount+1
	!print *, 'DRIVE_SEM::check5a'
            call savefield (Tdomain,ntime,rg,icount) 
	!print *, 'DRIVE_SEM::check5b'
	
         endif
      !endif
   enddo

endif

write (*,*) "Simulation is finished ",rg

!stop
write (*,*) "Deallocate fields ",rg
call deallocate_domain(Tdomain,rg)
!endif
T_comp2=MPI_WTIME()
if (rg==0) then
open (33,file="computation.time",form="formatted",status="unknown")
write(33,*) 'The total computation time is ', T_comp2-T_comp0, 'seconds'
write(33,*) 'in which ', T_comp1-T_comp0, 'seconds are spent for the field generation'
close (33)
endif
! Work in progress...
write (*,*) "END "


call mpi_barrier(MPI_COMM_WORLD,code)


call mpi_finalize(code)

end program drive_sem

 
! $Id: drive_sem.f90,v 1.3 2007/08/30 15:59:34 taquanganh Exp $ !
