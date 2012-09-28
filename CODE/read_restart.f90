subroutine read_restart (Tdomain,rg)

! Modified by Elise Delavaud (20/03/06)


use sdomain

implicit none

type (domain), intent (INOUT):: Tdomain
integer, intent (IN) :: rg

! local variables
integer :: n,ngllx,nglly,ngllz,i,j,k,ngll,ngll1,ngll2,ndt2
character (len=20) :: fnamef
logical :: finished_save,existe,StoreTrace

ndt2=75

 inquire (file="SBackup/status",EXIST=existe)
   if ( existe ) then 
      open (611,file="SBackup/status", status="old")
      read (611,*) finished_save
      close (611)
      if (finished_save) then
         write (fnamef,"(a,I2.2)") "SBackup/backup",rg
      else
         print *, 'NO BACKUP FOUND, compute from scratch'
         return
      endif
   else
      print *, 'NO BACKUP FOUND, compute from scratch'
      return
   endif

open (61,file=fnamef,status="unknown",form="formatted")
read(61,*) Tdomain%TimeD%rtime, Tdomain%TimeD%NtimeMin
print *,'read_restart::check0',rg
! Save Fields for Elements
do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
   if ( .not. Tdomain%specel(n)%PML ) then   
      do k = 1,ngllz-2
         do j = 1,nglly-2
            do i = 1,ngllx-2
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,0)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,1)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,2)
!!$               !U_{n} for TimeReverse
!!$               read(61,*) Tdomain%specel(n)%D0(i,j,k,0)
!!$               read(61,*) Tdomain%specel(n)%D0(i,j,k,1)
!!$               read(61,*) Tdomain%specel(n)%D0(i,j,k,2)
!!$
               !if (Tdomain%TimeReverse) Tdomain%specel(n)%Displ(i,j,k,:)=2*Tdomain%specel(n)%Displ(i,j,k,:)-Tdomain%specel(n)%D0(i,j,k,:)
!!$                if (Tdomain%TimeReverse) Tdomain%specel(n)%Veloc(i,j,k,:)=-Tdomain%specel(n)%Veloc(i,j,k,:)
            enddo
         enddo
      enddo
   else
      do k = 1,ngllz-2
       do j = 1,nglly-2
        do i = 1,ngllx-2
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,2)
        enddo
      enddo 
     enddo 
      do k = 0,ngllz-1
       do j = 0,nglly-1
        do i = 0,ngllx-1
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,2)

if (Tdomain%specel(n)%aniso) then
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress4b(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress5b(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Diagonal_Stress6b(i,j,k,2)

                  read(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress3(i,j,k,2)

                  read(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress4a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress4b(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress5b(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6a(i,j,k,2)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,0)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,1)
                  read(61,*) Tdomain%specel(n)%Residual_Stress6b(i,j,k,2)

               endif



        enddo
      enddo 
     enddo 
   endif 
enddo

print *,'read_restart::check1',rg

! Save Fields for Faces
do n = 0,Tdomain%n_face-1
   ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
   if (.not. Tdomain%sFace(n)%PML ) then
    do j = 1,ngll2-2
     do i = 1,ngll1-2
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,0)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,1)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,2)

!!$       !U_{n} for TimeReverse
!!$       read(61,*) Tdomain%sFace(n)%D0(i,j,0)
!!$       read(61,*) Tdomain%sFace(n)%D0(i,j,1)
!!$       read(61,*) Tdomain%sFace(n)%D0(i,j,2)
!!$       !if (Tdomain%TimeReverse) Tdomain%sFace(n)%Displ(i,j,:)=2*Tdomain%sFace(n)%Displ(i,j,:)-Tdomain%sFace(n)%D0(i,j,:)
!!$       if (Tdomain%TimeReverse) Tdomain%sFace(n)%Veloc(i,j,:)=-Tdomain%sFace(n)%Veloc(i,j,:)
    enddo
 enddo
else
 do j = 1,ngll2-2
      do i = 1,ngll1-2
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,2)
      enddo
     enddo
   endif
enddo
print *,'read_restart::check2',rg
! Save Fields for Edges
do n = 0,Tdomain%n_edge-1
   ngll = Tdomain%sEdge(n)%ngll
   if (.not. Tdomain%sEdge(n)%PML ) then
      do i = 1,ngll-2
         read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
         read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
         read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
         read(61,*) Tdomain%sEdge(n)%Displ(i,0)
         read(61,*) Tdomain%sEdge(n)%Displ(i,1)
         read(61,*) Tdomain%sEdge(n)%Displ(i,2)

!!$         !U_{n} for TimeReverse
!!$         read(61,*) Tdomain%sEdge(n)%D0(i,0)
!!$         read(61,*) Tdomain%sEdge(n)%D0(i,1)
!!$         read(61,*) Tdomain%sEdge(n)%D0(i,2)
!!$         !if (Tdomain%TimeReverse) Tdomain%sEdge(n)%Displ(i,:)=2*Tdomain%sEdge(n)%Displ(i,:)-Tdomain%sEdge(n)%D0(i,:)
!!$         if (Tdomain%TimeReverse) Tdomain%sEdge(n)%Veloc(i,:)=-Tdomain%sEdge(n)%Veloc(i,:)


      enddo
   else
     do i = 1,ngll-2
       read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,2)
     enddo
   endif
enddo
print *,'read_restart::check3',rg
! Save Fields for Vertices
do n = 0,Tdomain%n_vertex-1
   if (.not. Tdomain%sVertex(n)%PML ) then
      do i = 0,2
         read(61,*) Tdomain%sVertex(n)%Veloc(i)
         read(61,*) Tdomain%sVertex(n)%Displ(i)

!!$         !U_{n} for TimeReverse 
!!$         read(61,*) Tdomain%sVertex(n)%D0(i)
      enddo
!!$      if (Tdomain%TimeReverse) Tdomain%sVertex(n)%Displ(:)=2*Tdomain%sVertex(n)%Displ(:)-Tdomain%sVertex(n)%D0(:)
!!$ if (Tdomain%TimeReverse) Tdomain%sVertex(n)%Veloc(:)=- Tdomain%sVertex(n)%Veloc(:)
   else
      do i = 0,2
         read(61,*) Tdomain%sVertex(n)%Veloc(i)
         read(61,*) Tdomain%sVertex(n)%Veloc1(i)
         read(61,*) Tdomain%sVertex(n)%Veloc2(i)
         read(61,*) Tdomain%sVertex(n)%Veloc3(i)
      enddo
   endif
enddo

print *,'read_restart::check4',rg
!!$j=0
!Load StoreTrace
if (Tdomain%logicD%save_trace) then
   read (61,*)
   do n = 0, Tdomain%n_receivers-1 
      if (rg == Tdomain%sReceiver(n)%proc) then
         if ( Tdomain%sReceiver(n)%flag == 1)  then
            allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%TimeD%ntrace-1, 0:2) )
         elseif ( Tdomain%sReceiver(n)%flag == 2) then
            allocate (Tdomain%sReceiver(n)%StoreTrace (0:(Tdomain%TimeD%ntrace-1)/75, 0:2) )
         endif
         read (61,*) StoreTrace
         if (StoreTrace) then
            do i=0, size(Tdomain%sReceiver(n)%StoreTrace,1)-1
               read (61,*) Tdomain%sReceiver(n)%StoreTrace(i,0)
!!$j=j+1
!!$print *,j
               read (61,*) Tdomain%sReceiver(n)%StoreTrace(i,1)
!!$j=j+1
!!$print *,j
               read (61,*) Tdomain%sReceiver(n)%StoreTrace(i,2)
!!$j=j+1
!!$print *,j
            enddo
         endif
      endif
   enddo
endif


close(61) 



return
end subroutine read_restart
