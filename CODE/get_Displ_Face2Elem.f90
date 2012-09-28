subroutine get_Displ_Face2Elem (Tdomain,n)

! Modified by Elise Delavaud 26/01/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: nf,i,j,ngllx,nglly,ngllz,ngll1,ngll2
 
!print*,'hello'
ngllx = Tdomain%specel(n)%ngllx; nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

nf = Tdomain%specel(n)%near_faces(0)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(0) == 0 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb00!!!!!!!!!!'
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,0,0) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0) 
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,0,1) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,1) 
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,0,2) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(0) == 1 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb01!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,j,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,j,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,j,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 2 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb02!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(i,nglly-1-j,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(i,nglly-1-j,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(i,nglly-1-j,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 3 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb03!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 4 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb04!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,i,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,i,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,i,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 5 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb05!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,i,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,i,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,i,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 6 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb06!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,nglly-1-i,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,nglly-1-i,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,nglly-1-i,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 7 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb07!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF0'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif

nf = Tdomain%specel(n)%near_faces(1)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(1) == 0 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb10!!!!!!!!!!'
  Tdomain%specel(n)%Forces(1:ngllx-2,0,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0) 
  Tdomain%specel(n)%Forces(1:ngllx-2,0,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,1) 
  Tdomain%specel(n)%Forces(1:ngllx-2,0,1:ngllz-2,2) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(1) == 1 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb11!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,0,j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 2 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb12!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(i,0,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(i,0,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(i,0,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 3 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb13!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 4 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb14!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,0,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,0,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,0,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 5 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb15!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,0,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,0,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,0,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 6 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb16!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,0,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,0,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,0,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)  
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 7 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb17!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF1'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif

!if (n==9) then
!write(11,*) 'Elem9getF1'
!do i=0,4
!  do j=0,4
!    write(11,*) i,j,Tdomain%specel(9)%Forces(i,j,0,0)
!  enddo
!enddo
!endif

nf = Tdomain%specel(n)%near_faces(2)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(2) == 0 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb20!!!!!!!!!!'
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0) 
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,1) 
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,1:ngllz-2,2) =  Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(2) == 1 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb21!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 2 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb22!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 3 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb23!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 4 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb24!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,j,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,j,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,j,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 5 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb25!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 6 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb26!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 7 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb27!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF2'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif
!if (n==9) then
!write(11,*) 'Elem9getF2'
!do i=0,4
!  do j=0,4
!    write(11,*) i,j,Tdomain%specel(9)%Forces(i,j,0,0)
!  enddo
!enddo
!endif

nf = Tdomain%specel(n)%near_faces(3)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(3) == 0 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb30!!!!!!!!!!'
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:ngllz-2,0) 
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:ngllz-2,1) 
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,1:ngllz-2,2) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:ngllz-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(3) == 1 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb31!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 2 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb32!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 3 ) then
!if (.not. (ngll1==ngllx .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb33!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 4 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb34!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,nglly-1,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,nglly-1,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,nglly-1,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 5 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb35!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 6 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb36!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 7 ) then
!if (.not. (ngll1==ngllz .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb37!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF3'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif
!if (n==9) then
!write(11,*) 'Elem9getF3'
!do i=0,4
!  do j=0,4
!    write(11,*) i,j,Tdomain%specel(9)%Forces(i,j,0,0)
!  enddo
!enddo
!endif

nf = Tdomain%specel(n)%near_faces(4)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(4) == 0 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb40!!!!!!!!!!'
  Tdomain%specel(n)%Forces(0,1:nglly-2,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(1:nglly-2,1:ngllz-2,0) 
  Tdomain%specel(n)%Forces(0,1:nglly-2,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(1:nglly-2,1:ngllz-2,1) 
  Tdomain%specel(n)%Forces(0,1:nglly-2,1:ngllz-2,2) =  Tdomain%sFace(nf)%Forces(1:nglly-2,1:ngllz-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(4) == 1 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb41!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,nglly-1-i,j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,nglly-1-i,j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,nglly-1-i,j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 2 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb42!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 3 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllz) ) print*,'!!!!!!!!!!Pb43!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 4 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb44!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,j,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,j,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,j,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 5 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb45!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Forces(0,nglly-1-j,i,0) =  Tdomain%sFace(nf)%Forces(i,j,0)
    Tdomain%specel(n)%Forces(0,nglly-1-j,i,1) =  Tdomain%sFace(nf)%Forces(i,j,1)
    Tdomain%specel(n)%Forces(0,nglly-1-j,i,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 6 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb46!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 7 ) then
!if (.not. (ngll1==ngllz .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb47!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Forces(i,j,2) 
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF4'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif
!if (n==9) then
!write(11,*) 'Elem9getF4'
!do i=0,4
!  do j=0,4
!    write(11,*) i,j,Tdomain%specel(9)%Forces(i,j,0,0)
!  enddo
!enddo
!endif

nf = Tdomain%specel(n)%near_faces(5)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(5) == 0 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb50!!!!!!!!!!'
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:nglly-2,0) 
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:nglly-2,1) 
  Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(1:ngllx-2,1:nglly-2,2) 
else if ( Tdomain%specel(n)%Orient_Faces(5) == 1 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb51!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 2 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb52!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 3 ) then
!if (.not. (ngll1==ngllx .and. ngll2==nglly) ) print*,'!!!!!!!!!!Pb53!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 4 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb54!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,i,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,i,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,i,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 5 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb55!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 6 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb56!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 7 ) then
!if (.not. (ngll1==nglly .and. ngll2==ngllx) ) print*,'!!!!!!!!!!Pb57!!!!!!!!!!'
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(i,j,0) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(i,j,1) 
     Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,2) =  Tdomain%sFace(nf)%Forces(i,j,2)
   enddo
 enddo
endif

!if (n==13) then
!write(12,*) 'Elem13getF5'
!do i=0,4
!  do j=0,4
!    write(12,*) i,j,Tdomain%specel(13)%Forces(0,i,j,2)
!  enddo
!enddo
!endif
!if (n==9) then
!write(11,*) 'Elem9getF5'
!do i=0,4
!  do j=0,4
!    write(11,*) i,j,Tdomain%specel(9)%Forces(i,j,0,0)
!  enddo
!enddo
!endif

return
end subroutine get_Displ_Face2Elem
