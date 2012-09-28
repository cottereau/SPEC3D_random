subroutine save_checkpoint (Tdomain,rtime,it,rg)

! Modified by Elise Delavaud (20/03/06)


use MPI
use sdomain

implicit none

type (domain), intent (IN):: Tdomain
integer, intent (IN) :: it,rg
real, intent (IN) :: rtime

! local variables
integer :: n,ngllx,nglly,ngllz,i,j,k,ngll,ngll1,ngll2,code
character (len=100) :: fnamef,commande

   write (fnamef,"(a,I2.2)") "SBackup/backup",rg
   open (611,file="SBackup/status",status="replace",form="formatted")
   write(611,*),'.false.'
   close(611)

open (61,file=fnamef,status="unknown",form="formatted")
write(61,*) rtime,it

! Save Fields for Elements
do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
   if ( .not. Tdomain%specel(n)%PML ) then   
      do k = 1,ngllz-2
         do j = 1,nglly-2
            do i = 1,ngllx-2
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Displ(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Displ(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Displ(i,j,k,2)
  
!!$               ! U_{n} for TimeReverse
!!$               write(61,*) Tdomain%specel(n)%D0(i,j,k,0)
!!$               write(61,*) Tdomain%specel(n)%D0(i,j,k,1)
!!$               write(61,*) Tdomain%specel(n)%D0(i,j,k,2)
!!$               
            enddo
         enddo
      enddo
   else
      do k = 1,ngllz-2
         do j = 1,nglly-2
            do i = 1,ngllx-2
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,2)
            enddo
         enddo
      enddo
      do k = 0,ngllz-1
         do j = 0,nglly-1
            do i = 0,ngllx-1
               write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
               write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
               write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
               write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,2)
               if (Tdomain%specel(n)%aniso) then
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,2)

                  write(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,2)

                  write(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,2)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,0)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,1)
                  write(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,2)

               endif
            enddo
         enddo
      enddo
   endif
enddo
!print *, 'save_checkpoint::check1',rg
! Save Fields for Faces
do n = 0,Tdomain%n_face-1
   ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
   if (.not. Tdomain%sFace(n)%PML ) then
    do j = 1,ngll2-2
     do i = 1,ngll1-2
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       write(61,*) Tdomain%sFace(n)%Displ(i,j,0)
       write(61,*) Tdomain%sFace(n)%Displ(i,j,1)
       write(61,*) Tdomain%sFace(n)%Displ(i,j,2)
       
!!$       !U_{n} for TimeReverse
!!$       write(61,*) Tdomain%sFace(n)%D0(i,j,0)
!!$       write(61,*) Tdomain%sFace(n)%D0(i,j,1)
!!$       write(61,*) Tdomain%sFace(n)%D0(i,j,2)
!!$

      enddo
     enddo
   else
     do j = 1,ngll2-2
      do i = 1,ngll1-2
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       write(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       write(61,*) Tdomain%sFace(n)%Veloc1(i,j,0)
       write(61,*) Tdomain%sFace(n)%Veloc1(i,j,1)
       write(61,*) Tdomain%sFace(n)%Veloc1(i,j,2)
       write(61,*) Tdomain%sFace(n)%Veloc2(i,j,0)
       write(61,*) Tdomain%sFace(n)%Veloc2(i,j,1)
       write(61,*) Tdomain%sFace(n)%Veloc2(i,j,2)
       write(61,*) Tdomain%sFace(n)%Veloc3(i,j,0)
       write(61,*) Tdomain%sFace(n)%Veloc3(i,j,1)
       write(61,*) Tdomain%sFace(n)%Veloc3(i,j,2)
      enddo
     enddo
   endif
enddo
!print *, 'save_checkpoint::check2',rg
! Save Fields for Edges
do n = 0,Tdomain%n_edge-1
   ngll = Tdomain%sEdge(n)%ngll
   if (.not. Tdomain%sEdge(n)%PML ) then
     do i = 1,ngll-2
       write(61,*) Tdomain%sEdge(n)%Veloc(i,0)
       write(61,*) Tdomain%sEdge(n)%Veloc(i,1)
       write(61,*) Tdomain%sEdge(n)%Veloc(i,2)
       write(61,*) Tdomain%sEdge(n)%Displ(i,0)
       write(61,*) Tdomain%sEdge(n)%Displ(i,1)
       write(61,*) Tdomain%sEdge(n)%Displ(i,2)

!!$       !U_{n} for TimeReverse
!!$       write(61,*) Tdomain%sEdge(n)%D0(i,0)
!!$       write(61,*) Tdomain%sEdge(n)%D0(i,1)
!!$       write(61,*) Tdomain%sEdge(n)%D0(i,2)

     enddo
   else
     do i = 1,ngll-2
       write(61,*) Tdomain%sEdge(n)%Veloc(i,0)
       write(61,*) Tdomain%sEdge(n)%Veloc(i,1)
       write(61,*) Tdomain%sEdge(n)%Veloc(i,2)
       write(61,*) Tdomain%sEdge(n)%Veloc1(i,0)
       write(61,*) Tdomain%sEdge(n)%Veloc1(i,1)
       write(61,*) Tdomain%sEdge(n)%Veloc1(i,2)
       write(61,*) Tdomain%sEdge(n)%Veloc2(i,0)
       write(61,*) Tdomain%sEdge(n)%Veloc2(i,1)
       write(61,*) Tdomain%sEdge(n)%Veloc2(i,2)
       write(61,*) Tdomain%sEdge(n)%Veloc3(i,0)
       write(61,*) Tdomain%sEdge(n)%Veloc3(i,1)
       write(61,*) Tdomain%sEdge(n)%Veloc3(i,2)
     enddo
   endif
enddo
!print *, 'save_checkpoint::check3',rg
! Save Fields for Vertices
do n = 0,Tdomain%n_vertex-1
   if (.not. Tdomain%sVertex(n)%PML ) then
     do i = 0,2
       write(61,*) Tdomain%sVertex(n)%Veloc(i)
       write(61,*) Tdomain%sVertex(n)%Displ(i)

!!$       write(61,*) Tdomain%sVertex(n)%D0(i)
     enddo
   else
     do i = 0,2
       write(61,*) Tdomain%sVertex(n)%Veloc(i)
       write(61,*) Tdomain%sVertex(n)%Veloc1(i)
       write(61,*) Tdomain%sVertex(n)%Veloc2(i)
       write(61,*) Tdomain%sVertex(n)%Veloc3(i)
     enddo
   endif
enddo

!!$commande= 'wc -l ' // fnamef
!!$call system(commande)

!print *, 'save_checkpoint::check4',rg

!Save StoreTrace
if (Tdomain%logicD%save_trace) then
   write (61,*) "SAVE_STORETRACE"
   do n = 0, Tdomain%n_receivers-1 
      if (rg == Tdomain%sReceiver(n)%proc) then
         write (61,*) allocated(Tdomain%sReceiver(n)%StoreTrace)
         if (allocated(Tdomain%sReceiver(n)%StoreTrace)) then
            do i=0, size(Tdomain%sReceiver(n)%StoreTrace,1)-1
               write (61,*) Tdomain%sReceiver(n)%StoreTrace(i,0)
               write (61,*) Tdomain%sReceiver(n)%StoreTrace(i,1)
               write (61,*) Tdomain%sReceiver(n)%StoreTrace(i,2)
!!$if (rg==6) print *, Tdomain%sReceiver(n)%StoreTrace(i,0), Tdomain%sReceiver(n)%StoreTrace(i,1), Tdomain%sReceiver(n)%StoreTrace(i,2)
            enddo
         endif
      endif
   enddo
endif

close(61) 

!!$commande= 'wc -l ' // fnamef
!!$call system(commande)

call mpi_barrier(MPI_COMM_WORLD,code)

if (rg==0) then

   open (611,file="SBackup/status",status="replace",form="formatted")
   write(611,*),'.true.'
   close(611)

endif



return
end subroutine save_checkpoint
