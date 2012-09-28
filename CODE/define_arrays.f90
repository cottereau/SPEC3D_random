subroutine Define_Arrays(Tdomain, rg)

  ! Modified by Gaetano Festa 23/02/2005
  ! Modified by Paul Cupillard 12/11/2005
  ! Modified by Elise Delavaud 15/02/2006


  use sdomain
  use pig


  implicit none

  include 'mpif.h'

  type (domain), intent (INOUT), target :: Tdomain
  integer, intent(IN) :: rg

  integer :: n, mat, ngllx,nglly,ngllz, ngll1,ngll2, ngll,ngllPML,ngllSO, i,j,k, n_elem, nf,ne,nv,nv_aus, idef, code,notOpen
  integer :: shift, I_give_to, I_take_from, n_rings
  integer, parameter :: etiquette=100
  integer, dimension(mpi_status_size) :: statut
  real :: vp, ri,rj,rk, dx
  real, external :: pow
  real, dimension (:,:,:), allocatable :: xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz, Jac
  real, dimension (:,:,:), allocatable :: Rlam,Rmu,RKmod, Whei, LocMassMat, wx,wy,wz, Id, Store_Btn

  integer :: desc, SZG, nbprocs
  integer, dimension(:), allocatable :: NbElem, NbGlob

  real,dimension(:,:,:,:,:),allocatable, target::c
logical::existe,comparateur,from_scratch_format
integer::n_elem_proc,n_glob_proc
real :: maxRate


 notOpen=1
  !print *,rg,Tdomain%n_glob_points

!!$  call MPI_FILE_OPEN(MPI_COMM_WORLD,"RandomField.LAM",&
!!$       MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)
SZG=3*Tdomain%n_glob_points
call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
allocate(NbElem(0:nbprocs-1))
call MPI_ALLGATHER(Tdomain%n_elem,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
     MPI_COMM_WORLD,code)
!print *,rg,Tdomain%n_elem
allocate(NbGlob(0:nbprocs-1))
call MPI_ALLGATHER(Tdomain%n_glob_points,1,MPI_INTEGER,NbGlob,1,MPI_INTEGER, &
     MPI_COMM_WORLD,code)
!desc=-1
inquire (file="../RandomField.LAM",EXIST=existe)
1000 if (.not. (existe)) then
   if (Tdomain%SaveMecaField) then
      call MPI_FILE_OPEN(MPI_COMM_WORLD,"../RandomField.LAM",&
           MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,notOpen)
      call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_elem,1,MPI_INTEGER,statut,code)
      call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_glob_points,1,MPI_INTEGER,statut,code)
      call MPI_FILE_WRITE_ORDERED(desc,Tdomain%GlobCoord,SZG,MPI_DOUBLE_PRECISION,statut,code)
   endif
   from_scratch_format=.true.
else 
!print *, 'check_RANDOM.LAM'
   call MPI_FILE_OPEN(MPI_COMM_WORLD,"../RandomField.LAM",&
        MPI_MODE_RDWR,MPI_INFO_NULL,desc,code)
   call MPI_FILE_READ_ORDERED(desc,n_elem_proc,1,MPI_INTEGER,statut,code)
   call MPI_FILE_READ_ORDERED(desc,n_glob_proc,1,MPI_INTEGER,statut,code)
   comparateur=((n_elem_proc==Tdomain%n_elem) .and. (n_glob_proc==Tdomain%n_glob_points))
   call MPI_ALLREDUCE(comparateur,existe,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,code) 
   if (.not. (existe)) goto 1000
endif

!call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
!allocate(NbElem(0:nbprocs-1))
!call MPI_ALLGATHER(Tdomain%n_elem,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
!     MPI_COMM_WORLD,code)
!print *,rg,Tdomain%n_elem
!allocate(NbGlob(0:nbprocs-1))
!call MPI_ALLGATHER(Tdomain%n_glob_points,1,MPI_INTEGER,NbGlob,1,MPI_INTEGER, &
!     MPI_COMM_WORLD,code)

!print *, 'check3'

!from_scratch_format=.true.

if (Tdomain%SaveMecaField) then
   if (from_scratch_format) then
      do  n = 0,Tdomain%n_elem-1
         mat = Tdomain%specel(n)%mat_index
         ngllx = Tdomain%specel(n)%ngllx
         nglly = Tdomain%specel(n)%nglly
         ngllz = Tdomain%specel(n)%ngllz

         if (.not. allocated(c))  allocate(c(1:6,1:6,0:ngllx-1,0:nglly-1,0:ngllz-1))
         call define_mechanical_fields_aniso(Tdomain%sSubdomain(mat), &
              ngllx,nglly,ngllz,Tdomain%specel(n)%Iglobnum, &
              Tdomain%GlobCoord,c,Tdomain%specel(n)%PML, &
              desc,NbElem,NbGlob,nbprocs,n,rg,from_scratch_format,Tdomain%SaveMecaField)
      enddo
      from_scratch_format=.false.

      call MPI_BARRIER (MPI_COMM_WORLD, code) 
      call MPI_FILE_CLOSE(desc,code) 
      print *, 'define_array::The size of RandomField.LAM is reached'
      call system('ls -lstrh ../RandomField.LAM')
      existe=.true.
      deallocate(c)
      goto 1000
   endif
else
   from_scratch_format=.false.
endif

!print *, 'define_array',desc

!!! Attribute elastic properties from material !!!
  do  n = 0,Tdomain%n_elem-1
     !do n=125,125
     mat = Tdomain%specel(n)%mat_index
!!$     Tdomain%specel(n)%aniso = .false. 
!!$     if (Tdomain%sSubDomain(mat)%material_type(2:2) == "A") Tdomain%specel(n)%aniso = .true.


     ngllx = Tdomain%specel(n)%ngllx
     nglly = Tdomain%specel(n)%nglly
     ngllz = Tdomain%specel(n)%ngllz

     Tdomain%specel(n)%Density = Tdomain%sSubDomain(mat)%Ddensity
     !     Tdomain%specel(n)%Lambda = Tdomain%sSubDomain(mat)%DLambda
     !     Tdomain%specel(n)%Mu = Tdomain%sSubDomain(mat)%DMu

     allocate (Jac(0:ngllx-1,0:nglly-1,0:ngllz-1))

     allocate (xix (0:ngllx-1,0:nglly-1,0:ngllz-1)) 
     allocate (xiy (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (xiz (0:ngllx-1,0:nglly-1,0:ngllz-1))    

     allocate (etax (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (etay (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (etaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

     allocate (zetax (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (zetay (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (zetaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

     allocate (Whei (0:ngllx-1,0:nglly-1,0:ngllz-1))

     do k = 0, ngllz -1 
        do j = 0,nglly-1 
           do i = 0,ngllx-1
              Whei (i,j,k) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j) &
                   * Tdomain%sSubdomain(mat)%GLLwz(k)
           enddo
        enddo
     enddo
!!$print *,'define_array,xix000',xix(0,0,0)
!!$print *,'define array,Invgrad00000', Tdomain%specel(n)%InvGrad(0,0,0,0,0)
     xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)
     xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
     xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)

     etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
     etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
     etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

     zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
     zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
     zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)

     Jac  = Tdomain%specel(n)%Jacob
     if (Tdomain%specel(n)%aniso) then
        allocate(c(1:6,1:6,0:ngllx-1,0:nglly-1,0:ngllz-1))
	if (.not. (existe)) Tdomain%sSubdomain(mat)%reload_mat=.false.
!print *,'define_check0'
        call define_mechanical_fields_aniso(Tdomain%sSubdomain(mat), &
             ngllx,nglly,ngllz,Tdomain%specel(n)%Iglobnum, &
             Tdomain%GlobCoord,c,Tdomain%specel(n)%PML, &
             desc,NbElem,NbGlob,nbprocs,n,rg,from_scratch_format,Tdomain%SaveMecaField)
!print *,'define_check1'
!       if ((Tdomain%logicD%save_snapshots==.true.) .and.  (Tdomain%Field_Order(2)==.true.)) Tdomain%specel(n)%c=c
        if (((Tdomain%logicD%save_snapshots==.true.) .and.  (Tdomain%Field_Order(2)==.true.))  .or. &
             (Tdomain%logicD%save_energy==.true.)) Tdomain%specel(n)%c=>c
!print *, 'check1_Tdomain%specel(n)%c(6,6,1,1,1)%',Tdomain%specel(n)%c(6,6,1,1,1)
     else
        allocate (RKmod (0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate (Rlam(0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate (Rmu(0:ngllx-1,0:nglly-1,0:ngllz-1))
        call define_mechanical_fields_iso(Tdomain%sSubdomain(mat), &
             ngllx,nglly,ngllz,Tdomain%specel(n)%Iglobnum, &
             Tdomain%GlobCoord,RLam,RMu,Tdomain%specel(n)%PML, &
             desc,NbElem,NbGlob,nbprocs,n,rg)
        Tdomain%specel(n)%Lambda = RLam
        Tdomain%specel(n)%Mu = RMu
        RKmod = Rlam + 2. * Rmu
        !if ((Tdomain%logicD%save_snapshots==.true.) .and.  (Tdomain%Field_Order(2)==.true.)) then
         !  Tdomain%specel(n)%c=0.
         !  Tdomain%specel(n)%c(1:3,1:3,:,:,:)=RLam
         !  Tdomain%specel(n)%c(1,1,:,:,:)=2*RMu+Tdomain%specel(n)%c(1,1,:,:,:)
         !  Tdomain%specel(n)%c(2,2,:,:,:)=2*RMu+Tdomain%specel(n)%c(2,2,:,:,:)
         !  Tdomain%specel(n)%c(3,3,:,:,:)=2*RMu+Tdomain%specel(n)%c(3,3,:,:,:)
         !  Tdomain%specel(n)%c(4,4,:,:,:)=RMu
         !  Tdomain%specel(n)%c(5,5,:,:,:)=RMu
         ! Tdomain%specel(n)%c(6,6,:,:,:)=RMu
        !endif
     endif
     Tdomain%specel(n)%MassMat = Whei*Tdomain%specel(n)%Density*Jac
     ! Acoeff. Computation
     if (.not. Tdomain%specel(n)%aniso) then !Isotropy
        if (.not. Tdomain%specel(n)%PML) then ! None-PML isotropic

           Tdomain%specel(n)%Acoeff(:,:,:,0) = -Whei*(RKmod*xix**2+Rmu*(xiy**2+xiz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,1) = -Whei*(RKmod*xix*etax+Rmu*(xiy*etay+xiz*etaz))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,2) = -Whei*(RKmod*xix*zetax+Rmu*(xiy*zetay+xiz*zetaz))*Jac 
           Tdomain%specel(n)%Acoeff(:,:,:,3) = -Whei*(Rlam+Rmu)*xix*xiy*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,4) = -Whei*(Rlam*xix*etay+Rmu*xiy*etax)*Jac  
           Tdomain%specel(n)%Acoeff(:,:,:,5) = -Whei*(Rlam*xix*zetay+Rmu*xiy*zetax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,6) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,7) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,8) = -Whei*(Rlam*xix*zetaz+rmu*xiz*zetax)*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,9) = -Whei*(RKmod*etax**2+Rmu* (etay**2+etaz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,10) = -Whei*(RKmod*etax*zetax+Rmu* (etay*zetay+etaz*zetaz))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,11) = -Whei*(Rlam*etax*xiy+Rmu*etay*xix)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,12) = -Whei*(Rlam+Rmu)*etay*etax*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,13) = -Whei*(Rlam*etax*zetay+Rmu*etay*zetax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,14) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,15) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,16) = -Whei*(Rlam*etax*zetaz+Rmu*etaz*zetax)*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,17) = -Whei*(RKmod*zetax**2+Rmu*  (zetay**2+zetaz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,18) = -Whei*(Rlam*zetax*xiy+Rmu*zetay*xix)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,19) = -Whei*(Rlam*zetax*etay+Rmu*zetay*etax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,20) = -Whei*(Rlam+Rmu)*zetax*zetay*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,21) = -Whei*(Rlam*zetax*xiz+Rmu*zetaz*xix)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,22) = -Whei*(Rlam*zetax*etaz+Rmu*zetaz*etax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,23) = -Whei*(Rlam+Rmu)*zetax*zetaz*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,24) = -Whei*(RKmod*xiy**2+Rmu* (xix**2+xiz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,25) = -Whei*(RKmod*xiy*etay+Rmu* (xix*etax+xiz*etaz))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,26) = -Whei*(RKmod*xiy*zetay+Rmu*  (xix*zetax+xiz*zetaz))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei*(Rlam+Rmu)*xiy*xiz*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei*(Rlam*etaz*xiy+Rmu*etay*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei*(Rlam*zetaz*xiy+Rmu*zetay*xiz)*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei*(RKmod*etay**2+Rmu* (etax**2+etaz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei*(RKmod*zetay*etay+Rmu* (zetax*etax+zetaz*etaz))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei*(Rlam*etay*xiz+Rmu*etaz*xiy)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei*(Rlam+Rmu)*etay*etaz*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei*(Rlam*zetaz*etay+Rmu*zetay*etaz)*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei*(RKmod*zetay**2+Rmu* (zetax**2+zetaz**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,36) = -Whei*(Rlam*xiz*zetay+Rmu*xiy*zetaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,37) = -Whei*(Rlam*zetay*etaz+Rmu*zetaz*etay)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,38) = -Whei*(Rlam+Rmu)*zetay*zetaz*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,39) = -Whei*(RKmod*xiz**2+Rmu*  (xix**2+xiy**2))*Jac 
           Tdomain%specel(n)%Acoeff(:,:,:,40) = -Whei*(RKmod*xiz*etaz+Rmu* (xix*etax+xiy*etay))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,41) = -Whei*(RKmod*xiz*zetaz+Rmu* (xix*zetax+xiy*zetay))*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,42) = -Whei*(RKmod*etaz**2+Rmu* (etax**2+etay**2))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,43) = -Whei*(RKmod*zetaz*etaz+Rmu* (zetax*etax+zetay*etay))*Jac

           Tdomain%specel(n)%Acoeff(:,:,:,44) = -Whei*(RKmod*zetaz**2+Rmu* (zetax**2+zetay**2))*Jac

        else !PML isotropic

           Tdomain%specel(n)%Acoeff(:,:,:,0) = RKmod *xix
           Tdomain%specel(n)%Acoeff(:,:,:,1) = RKmod *etax
           Tdomain%specel(n)%Acoeff(:,:,:,2) = RKmod *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,3) = RLam *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,4) = RLam *etay
           Tdomain%specel(n)%Acoeff(:,:,:,5) = RLam *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,6) = RLam *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,7) = RLam *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,8) = RLam *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,9) = RLam *xix
           Tdomain%specel(n)%Acoeff(:,:,:,10) = RLam *etax
           Tdomain%specel(n)%Acoeff(:,:,:,11) = RLam *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,12) = RKmod *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,13) = RKmod *etay
           Tdomain%specel(n)%Acoeff(:,:,:,14) = RKmod *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,15) = RKmod *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,16) = RKmod *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,17) = RKmod *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,18) = RMu *xix
           Tdomain%specel(n)%Acoeff(:,:,:,19) = RMu *etax
           Tdomain%specel(n)%Acoeff(:,:,:,20) = RMu *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,21) = RMu *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,22) = RMu *etay
           Tdomain%specel(n)%Acoeff(:,:,:,23) = RMu *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,24) = RMu *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,25) = RMu *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,26) = RMu *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei * xix * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei * xiy * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei * xiz * Jac

           Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei * etax * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei * etay * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei * etaz * Jac

           Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei * zetax * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei * zetay * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei * zetaz * Jac
        endif
     else !Anisotropy
        if (.not. Tdomain%specel(n)%PML) then !None-PML anisotropic
!!$        N.B : L'odre etait (au debut de l'adaptation de l'anisotropie) pour "sigma(...,4,5,6)=sigma(...,23,13,12)" mais il est pour "sigma(...,12,13,23)" dans d'autres parties du code donc, il faudrait echanger le role des indices "4" et "6" (ca y est !)  
           Tdomain%specel(n)%Acoeff(:,:,:,0) = -Whei*(c(1,1,:,:,:)*xix**2 + c(4,4,:,:,:)*xiy**2 + c(5,5,:,:,:)*xiz**2 + 2*c(1,4,:,:,:)*xix*xiy + &
                2*c(1,5,:,:,:)*xix*xiz + 2*c(4,5,:,:,:)*xiy*xiz)*Jac 
           Tdomain%specel(n)%Acoeff(:,:,:,1) = -Whei*(c(1,1,:,:,:)*xix*etax + c(4,4,:,:,:)*xiy*etay + c(5,5,:,:,:)*xiz*etaz + c(1,4,:,:,:)*&
                (xix*etay+xiy*etax) + c(1,5,:,:,:)*(xix*etaz+xiz*etax) + c(4,5,:,:,:)*(xiy*etaz+xiz*etay))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,2) = -Whei*(c(1,1,:,:,:)*xix*zetax + c(4,4,:,:,:)*xiy*zetay + c(5,5,:,:,:)*xiz*zetaz + c(1,4,:,:,:)*&
                (xix*zetay+xiy*zetax) + c(1,5,:,:,:)*(xix*zetaz+xiz*zetax) + c(4,5,:,:,:)*(xiy*zetaz+xiz*zetay))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,3) = -Whei*(c(1,4,:,:,:)*xix**2 + c(2,4,:,:,:)*xiy**2 + c(5,6,:,:,:)*xiz**2 + (c(1,2,:,:,:)+c(4,4,:,:,:))&
                *xix*xiy + (c(1,6,:,:,:)+c(4,5,:,:,:))*xix*xiz + (c(4,6,:,:,:)+c(2,5,:,:,:))*xiy*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,4) = -Whei*(c(1,2,:,:,:)*xix*etay + c(1,6,:,:,:)*xix*etaz + c(1,4,:,:,:)*xix*etax + c(2,4,:,:,:)*xiy*etay + &
                c(4,6,:,:,:)*xiy*etaz + c(4,4,:,:,:)*xiy*etax + c(2,5,:,:,:)*xiz*etay +c(5,6,:,:,:)*xiz*etaz + c(4,5,:,:,:)*xiz*etax)*Jac  
           Tdomain%specel(n)%Acoeff(:,:,:,5) = -Whei*(c(1,2,:,:,:)*xix*zetay + c(1,6,:,:,:)*xix*zetaz + c(1,4,:,:,:)*xix*zetax + c(2,4,:,:,:)*xiy*zetay &
                + c(4,6,:,:,:)*xiy*zetaz + c(4,4,:,:,:)*xiy*zetax + c(2,5,:,:,:)*xiz*zetay + c(5,6,:,:,:)*xiz*zetaz + c(4,5,:,:,:)*xiz*zetax)*Jac  
           Tdomain%specel(n)%Acoeff(:,:,:,6) = -Whei*(c(1,5,:,:,:)*xix**2 + c(4,6,:,:,:)*xiy**2 + c(3,5,:,:,:)*xiz**2 + (c(1,6,:,:,:)+c(4,5,:,:,:))&
                *xix*xiy + (c(1,3,:,:,:)+c(5,5,:,:,:))*xix*xiz + (c(3,4,:,:,:)+c(5,6,:,:,:))*xiy*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,7) = -Whei*(c(1,3,:,:,:)*xix*etaz + c(1,6,:,:,:)*xix*etay + c(1,5,:,:,:)*xix*etax + c(3,4,:,:,:)*xiy*etaz&
                + c(4,6,:,:,:)*xiy*etay + c(4,5,:,:,:)*xiy*etax + c(3,5,:,:,:)*xiz*etaz + c(5,6,:,:,:)*xiz*etay + c(5,5,:,:,:)*xiz*etax)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,8) = -Whei*(c(1,3,:,:,:)*xix*zetaz + c(1,6,:,:,:)*xix*zetay + c(1,5,:,:,:)*xix*zetax + c(3,4,:,:,:)*xiy*zetaz&
                + c(4,6,:,:,:)*xiy*zetay + c(4,5,:,:,:)*xiy*zetax + c(3,5,:,:,:)*xiz*zetaz + c(5,6,:,:,:)*xiz*zetay + c(5,5,:,:,:)*xiz*zetax)*Jac 
           !-------------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,9) = -Whei*(c(1,1,:,:,:)*etax**2 + c(4,4,:,:,:)*etay**2 + c(5,5,:,:,:)*etaz**2 + 2*c(1,4,:,:,:)*etax*etay&
                + 2*c(1,5,:,:,:)*etax*etaz + 2*c(4,5,:,:,:)*etay*etaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,10) = -Whei*(c(1,1,:,:,:)*etax*zetax + c(4,4,:,:,:)*etay*zetay + c(5,5,:,:,:)*etaz*zetaz + c(1,5,:,:,:)*&
                (etax*zetaz+etaz*zetax) + c(1,4,:,:,:)*(etax*zetay+etay*zetax) + c(4,5,:,:,:)*(etay*zetaz+etaz*zetay) )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,11) = -Whei*(c(1,2,:,:,:)*xiy*etax + c(1,6,:,:,:)*xiz*etax + c(1,4,:,:,:)*xix*etax + c(2,4,:,:,:)*xiy*etay +&
                c(4,6,:,:,:)*xiz*etay + c(4,4,:,:,:)*xix*etay + c(2,5,:,:,:)*xiy*etaz +c(5,6,:,:,:)*xiz*etaz + c(4,5,:,:,:)*xix*etaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,12) = -Whei*(c(1,4,:,:,:)*etax**2 + c(2,4,:,:,:)*etay**2 + c(5,6,:,:,:)*etaz**2 + (c(1,2,:,:,:)+c(4,4,:,:,:))&
                *etax*etay + (c(2,5,:,:,:)+c(4,6,:,:,:))*etay*etaz  + (c(1,6,:,:,:)+c(4,5,:,:,:))*etax*etaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,13) = -Whei*(c(1,2,:,:,:)*zetay*etax + c(1,6,:,:,:)*zetaz*etax + c(1,4,:,:,:)*zetax*etax + c(2,4,:,:,:)&
                *zetay*etay + c(4,6,:,:,:)*zetaz*etay + c(4,4,:,:,:)*zetax*etay + c(2,5,:,:,:)*zetay*etaz +c(5,6,:,:,:)*zetaz*etaz + c(4,5,:,:,:)*zetax*etaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,14) = -Whei*(c(1,6,:,:,:)*xiy*etax + c(1,3,:,:,:)*xiz*etax + c(1,5,:,:,:)*xix*etax + c(4,6,:,:,:)*xiy*etay + &
                c(3,4,:,:,:)*xiz*etay +c(4,5,:,:,:)*xix*etay + c(5,6,:,:,:)*xiy*etaz +c(3,5,:,:,:)*xiz*etaz + c(5,5,:,:,:)*xix*etaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,15) = -Whei*((c(1,6,:,:,:)+c(4,5,:,:,:))*etay*etax + (c(1,3,:,:,:)+c(5,5,:,:,:))*etaz*etax + (c(3,4,:,:,:)+&
                c(5,6,:,:,:))*etaz*etay  + c(1,5,:,:,:)*etax**2 + c(4,6,:,:,:)*etay**2 + c(3,5,:,:,:)*etaz**2)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,16) = -Whei*(c(1,6,:,:,:)*zetay*etax + c(1,3,:,:,:)*zetaz*etax + c(1,5,:,:,:)*zetax*etax + c(4,6,:,:,:)&
                *zetay*etay + c(3,4,:,:,:)*zetaz*etay +c(4,5,:,:,:)*zetax*etay + c(5,6,:,:,:)*zetay*etaz +c(3,5,:,:,:)*zetaz*etaz + c(5,5,:,:,:)*zetax*etaz )*Jac
           !------------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,17) = -Whei*(c(1,1,:,:,:)*zetax**2 + c(4,4,:,:,:)*zetay**2 + c(5,5,:,:,:)*zetaz**2 + 2*c(1,4,:,:,:)*&
                zetax*zetay + 2*c(1,5,:,:,:)*zetax*zetaz + 2*c(4,5,:,:,:)*zetay*zetaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,18) = -Whei*(c(1,2,:,:,:)*xiy*zetax + c(1,6,:,:,:)*xiz*zetax + c(1,4,:,:,:)*xix*zetax + c(2,4,:,:,:)&
                *xiy*zetay + c(4,6,:,:,:)*xiz*zetay +c(4,4,:,:,:)*xix*zetay + c(2,5,:,:,:)*xiy*zetaz +c(5,6,:,:,:)*xiz*zetaz + c(4,5,:,:,:)*xix*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,19) = -Whei*(c(1,2,:,:,:)*etay*zetax + c(1,6,:,:,:)*etaz*zetax + c(1,4,:,:,:)*etax*zetax + c(2,4,:,:,:)*etay*&
                zetay + c(4,6,:,:,:)*etaz*zetay +c(4,4,:,:,:)*etax*zetay + c(2,5,:,:,:)*etay*zetaz +c(5,6,:,:,:)*etaz*zetaz + c(4,5,:,:,:)*etax*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,20) = -Whei*((c(1,2,:,:,:)+c(4,4,:,:,:))*zetax*zetay + (c(1,6,:,:,:)+ c(4,5,:,:,:))*zetax*zetaz + (c(2,5,:,:,:)+&
                c(4,6,:,:,:))*zetay*zetaz + c(1,4,:,:,:)*zetax**2 + c(2,4,:,:,:)*zetay**2 +c(5,6,:,:,:)*zetaz**2)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,21) = -Whei*(c(1,6,:,:,:)*xiy*zetax + c(1,3,:,:,:)*xiz*zetax + c(1,5,:,:,:)*xix*zetax + c(4,6,:,:,:)*xiy*zetay +&
                c(3,4,:,:,:)*xiz*zetay + c(4,5,:,:,:)*xix*zetay + c(5,6,:,:,:)*xiy*zetaz +c(3,5,:,:,:)*xiz*zetaz + c(5,5,:,:,:)*xix*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,22) = -Whei*(c(1,6,:,:,:)*etay*zetax + c(1,3,:,:,:)*etaz*zetax + c(1,5,:,:,:)*etax*zetax + c(4,6,:,:,:)*etay*&
                zetay + c(3,4,:,:,:)*etaz*zetay +c(4,5,:,:,:)*etax*zetay + c(5,6,:,:,:)*etay*zetaz +c(3,5,:,:,:)*etaz*zetaz + c(5,5,:,:,:)*etax*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,23) = -Whei*((c(1,6,:,:,:)+c(4,5,:,:,:))*zetax*zetay + (c(1,3,:,:,:)+c(5,5,:,:,:))*zetax*zetaz + (c(3,4,:,:,:)+&
                c(5,6,:,:,:))*zetay*zetaz + c(1,5,:,:,:)*zetax**2 + c(4,6,:,:,:)*zetay**2 + c(3,5,:,:,:)*zetaz**2)*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,24) = -Whei*(c(4,4,:,:,:)*xix**2 + c(2,2,:,:,:)*xiy**2 + c(6,6,:,:,:)*xiz**2 + 2*c(2,4,:,:,:)*xix*xiy + &
                2*c(4,6,:,:,:)*xix*xiz + 2*c(2,6,:,:,:)*xiy*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,25) = -Whei*(c(4,4,:,:,:)*xix*etax + c(2,2,:,:,:)*xiy*etay + c(6,6,:,:,:)*xiz*etaz + c(2,4,:,:,:)*&
                (xix*etay+xiy*etax) + c(4,6,:,:,:)*(xix*etaz+xiz*etax) + c(2,6,:,:,:)*(xiy*etaz+xiz*etay))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,26) = -Whei*(c(4,4,:,:,:)*xix*zetax + c(2,2,:,:,:)*xiy*zetay + c(6,6,:,:,:)*xiz*zetaz + c(2,4,:,:,:)*&
                (xix*zetay+xiy*zetax) + c(4,6,:,:,:)*(xix*zetaz+xiz*zetax) + c(2,6,:,:,:)*(xiy*zetaz+xiz*zetay))*Jac 
           Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei*(c(4,5,:,:,:)*xix**2 + c(2,6,:,:,:)*xiy**2 + c(3,6,:,:,:)*xiz**2 + (c(2,5,:,:,:)+c(4,6,:,:,:))*&
                xix*xiy + (c(5,6,:,:,:)+c(3,4,:,:,:))*xix*xiz + (c(6,6,:,:,:)+c(2,3,:,:,:))*xiy*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei*(c(4,6,:,:,:)*etay*xix + c(3,4,:,:,:)*etaz*xix + c(4,5,:,:,:)*etax*xix + c(2,6,:,:,:)*etay*xiy + &
                c(2,3,:,:,:)*etaz*xiy + c(2,5,:,:,:)*etax*xiy + c(6,6,:,:,:)*etay*xiz +c(3,6,:,:,:)*etaz*xiz + c(5,6,:,:,:)*etax*xiz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei*(c(4,6,:,:,:)*zetay*xix + c(3,4,:,:,:)*zetaz*xix + c(4,5,:,:,:)*zetax*xix + c(2,6,:,:,:)*zetay*xiy +&
                c(2,3,:,:,:)*zetaz*xiy + c(2,5,:,:,:)*zetax*xiy + c(6,6,:,:,:)*zetay*xiz +c(3,6,:,:,:)*zetaz*xiz + c(5,6,:,:,:)*zetax*xiz )*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei*(c(4,4,:,:,:)*etax**2 + c(2,2,:,:,:)*etay**2 + c(6,6,:,:,:)*etaz**2 + 2*c(2,4,:,:,:)*etax*etay + &
                2*c(4,6,:,:,:)*etax*etaz + 2*c(2,6,:,:,:)*etay*etaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei*(c(4,4,:,:,:)*etax*zetax + c(2,2,:,:,:)*etay*zetay + c(6,6,:,:,:)*etaz*zetaz + c(2,4,:,:,:)*&
                (etax*zetay+etay*zetax) + c(4,6,:,:,:)*(etax*zetaz+etaz*zetax) + c(2,6,:,:,:)*(etay*zetaz+etaz*zetay))*Jac 
           Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei*(c(4,6,:,:,:)*xiy*etax + c(3,4,:,:,:)*xiz*etax + c(4,5,:,:,:)*xix*etax + c(2,6,:,:,:)*xiy*etay + &
                c(2,3,:,:,:)*xiz*etay + c(2,5,:,:,:)*xix*etay + c(6,6,:,:,:)*xiy*etaz +c(3,6,:,:,:)*xiz*etaz + c(5,6,:,:,:)*xix*etaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei*(c(4,5,:,:,:)*etax**2 + c(2,6,:,:,:)*etay**2 + c(3,6,:,:,:)*etaz**2 + (c(2,5,:,:,:)+c(4,6,:,:,:))*&
                etax*etay + (c(5,6,:,:,:)+c(3,4,:,:,:))*etax*etaz + (c(6,6,:,:,:)+c(2,3,:,:,:))*etay*etaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei*(c(4,6,:,:,:)*zetay*etax + c(3,4,:,:,:)*zetaz*etax + c(4,5,:,:,:)*zetax*etax + c(2,6,:,:,:)* zetay*&
                etay + c(2,3,:,:,:)*zetaz*etay +  c(2,5,:,:,:)*zetax*etay + c(6,6,:,:,:)*zetay*etaz +c(3,6,:,:,:)*zetaz*etaz + c(5,6,:,:,:)*zetax*etaz )*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei*(c(4,4,:,:,:)*zetax**2 + c(2,2,:,:,:)*zetay**2 + c(6,6,:,:,:)*zetaz**2 + 2*c(2,4,:,:,:)*zetax*zetay+&
                2*c(4,6,:,:,:)*zetax*zetaz + 2*c(2,6,:,:,:)*zetay*zetaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,36) = -Whei*(c(4,6,:,:,:)*xiy*zetax + c(3,4,:,:,:)*xiz*zetax + c(4,5,:,:,:)*xix*zetax + c(2,6,:,:,:)*xiy*zetay + &
                c(2,3,:,:,:)*xiz*zetay +c(2,5,:,:,:)*xix*zetay + c(6,6,:,:,:)*xiy*zetaz +c(3,6,:,:,:)*xiz*zetaz + c(5,6,:,:,:)*xix*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,37) = -Whei*(c(4,6,:,:,:)*etay*zetax + c(3,4,:,:,:)*etaz*zetax + c(4,5,:,:,:)*etax*zetax + c(2,6,:,:,:)*etay*&
                zetay + c(2,3,:,:,:)*etaz*zetay +c(2,5,:,:,:)*etax*zetay + c(6,6,:,:,:)*etay*zetaz +c(3,6,:,:,:)*etaz*zetaz + c(5,6,:,:,:)*etax*zetaz )*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,38) = -Whei*(c(4,5,:,:,:)*zetax**2 + c(2,6,:,:,:)*zetay**2 + c(3,6,:,:,:)*zetaz**2 + (c(2,5,:,:,:)+c(4,6,:,:,:))*&
                zetax*zetay + (c(5,6,:,:,:)+c(3,4,:,:,:))*zetax*zetaz + (c(6,6,:,:,:)+c(2,3,:,:,:))*zetay*zetaz)*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,39) = -Whei*(c(5,5,:,:,:)*xix**2 + c(6,6,:,:,:)*xiy**2 + c(3,3,:,:,:)*xiz**2 + 2*c(5,6,:,:,:)*xix*xiy + &
                2*c(3,5,:,:,:)*xix*xiz + 2*c(3,6,:,:,:)*xiy*xiz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,40) = -Whei*(c(5,5,:,:,:)*xix*etax + c(6,6,:,:,:)*xiy*etay + c(3,3,:,:,:)*xiz*etaz + c(5,6,:,:,:)*(xix*etay+xiy*&
                etax) + c(3,5,:,:,:)*(xix*etaz+xiz*etax) + c(3,6,:,:,:)*(xiy*etaz+xiz*etay))*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,41) = -Whei*(c(5,5,:,:,:)*xix*zetax + c(6,6,:,:,:)*xiy*zetay + c(3,3,:,:,:)*xiz*zetaz + c(5,6,:,:,:)*&
                (xix*zetay+xiy*zetax) + c(3,5,:,:,:)*(xix*zetaz+xiz*zetax) + c(3,6,:,:,:)*(xiy*zetaz+xiz*zetay))*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,42) = -Whei*(c(5,5,:,:,:)*etax**2 + c(6,6,:,:,:)*etay**2 + c(3,3,:,:,:)*etaz**2 + 2*c(5,6,:,:,:)*etax*etay + &
                2*c(3,5,:,:,:)*etax*etaz + 2*c(3,6,:,:,:)*etay*etaz)*Jac
           Tdomain%specel(n)%Acoeff(:,:,:,43) = -Whei*(c(5,5,:,:,:)*etax*zetax + c(6,6,:,:,:)*etay*zetay + c(3,3,:,:,:)*etaz*zetaz + c(5,6,:,:,:)*&
                (etax*zetay+etay*zetax) + c(3,5,:,:,:)*(etax*zetaz+etaz*zetax) + c(3,6,:,:,:)*(etay*zetaz+etaz*zetay))*Jac
           !-----------------------------------------------------
           Tdomain%specel(n)%Acoeff(:,:,:,44) = -Whei*(c(5,5,:,:,:)*zetax**2 + c(6,6,:,:,:)*zetay**2 + c(3,3,:,:,:)*zetaz**2 + 2*c(5,6,:,:,:)*zetax*zetay +&
                2*c(3,5,:,:,:)*zetax*zetaz + 2*c(3,6,:,:,:)*zetay*zetaz)*Jac


        else !PML anisotropic
!!!QA print **, n, 'PML',rg
           Tdomain%specel(n)%Acoeff(:,:,:,0) = c(1,1,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,1) = c(1,1,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,2) = c(1,1,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,3) = c(1,2,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,4) = c(1,2,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,5) = c(1,2,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,6) = c(1,3,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,7) = c(1,3,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,8) = c(1,3,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,9) = c(1,4,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,10) = c(1,4,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,11) = c(1,4,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,12) = c(1,4,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,13) = c(1,4,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,14) = c(1,4,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,15) = c(1,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,16) = c(1,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,17) = c(1,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,18) = c(1,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,19) = c(1,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,20) = c(1,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,21) = c(1,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,22) = c(1,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,23) = c(1,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,24) = c(1,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,25) = c(1,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,26) = c(1,6,:,:,:) *zetay
           !---------------------------------------------------- 1
           Tdomain%specel(n)%Acoeff(:,:,:,27) = c(1,2,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,28) = c(1,2,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,29) = c(1,2,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,30) = c(2,2,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,31) = c(2,2,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,32) = c(2,2,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,33) = c(2,3,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,34) = c(2,3,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,35) = c(2,3,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,36) = c(2,4,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,37) = c(2,4,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,38) = c(2,4,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,39) = c(2,4,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,40) = c(2,4,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,41) = c(2,4,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,42) = c(2,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,43) = c(2,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,44) = c(2,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,45) = c(2,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,46) = c(2,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,47) = c(2,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,48) = c(2,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,49) = c(2,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,50) = c(2,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,51) = c(2,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,52) = c(2,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,53) = c(2,6,:,:,:) *zetay
           !------------------------------------------------- 2
           Tdomain%specel(n)%Acoeff(:,:,:,54) = c(1,3,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,55) = c(1,3,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,56) = c(1,3,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,57) = c(2,3,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,58) = c(2,3,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,59) = c(2,3,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,60) = c(3,3,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,61) = c(3,3,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,62) = c(3,3,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,63) = c(3,4,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,64) = c(3,4,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,65) = c(3,4,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,66) = c(3,4,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,67) = c(3,4,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,68) = c(3,4,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,69) = c(3,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,70) = c(3,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,71) = c(3,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,72) = c(3,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,73) = c(3,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,74) = c(3,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,75) = c(3,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,76) = c(3,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,77) = c(3,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,78) = c(3,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,79) = c(3,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,80) = c(3,6,:,:,:) *zetay
           !------------------------------------------------- 3
           Tdomain%specel(n)%Acoeff(:,:,:,81) = c(1,4,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,82) = c(1,4,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,83) = c(1,4,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,84) = c(2,4,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,85) = c(2,4,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,86) = c(2,4,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,87) = c(3,4,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,88) = c(3,4,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,89) = c(3,4,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,90) = c(4,4,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,91) = c(4,4,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,92) = c(4,4,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,93) = c(4,4,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,94) = c(4,4,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,95) = c(4,4,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,96) = c(4,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,97) = c(4,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,98) = c(4,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,99) = c(4,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,100) = c(4,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,101) = c(4,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,102) = c(4,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,103) = c(4,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,104) = c(4,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,105) = c(4,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,106) = c(4,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,107) = c(4,6,:,:,:) *zetay
           !--------------------------------------------------- 4
           Tdomain%specel(n)%Acoeff(:,:,:,108) = c(1,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,109) = c(1,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,110) = c(1,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,111) = c(2,5,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,112) = c(2,5,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,113) = c(2,5,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,114) = c(3,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,115) = c(3,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,116) = c(3,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,117) = c(4,5,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,118) = c(4,5,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,119) = c(4,5,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,120) = c(4,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,121) = c(4,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,122) = c(4,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,123) = c(5,5,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,124) = c(5,5,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,125) = c(5,5,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,126) = c(5,5,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,127) = c(5,5,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,128) = c(5,5,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,129) = c(5,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,130) = c(5,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,131) = c(5,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,132) = c(5,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,133) = c(5,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,134) = c(5,6,:,:,:) *zetay
           !-------------------------------------------------  5
           Tdomain%specel(n)%Acoeff(:,:,:,135) = c(1,6,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,136) = c(1,6,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,137) = c(1,6,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,138) = c(2,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,139) = c(2,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,140) = c(2,6,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,141) = c(3,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,142) = c(3,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,143) = c(3,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,144) = c(4,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,145) = c(4,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,146) = c(4,6,:,:,:) *zetay

           Tdomain%specel(n)%Acoeff(:,:,:,147) = c(4,6,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,148) = c(4,6,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,149) = c(4,6,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,150) = c(5,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,151) = c(5,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,152) = c(5,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,153) = c(5,6,:,:,:) *xix
           Tdomain%specel(n)%Acoeff(:,:,:,154) = c(5,6,:,:,:) *etax
           Tdomain%specel(n)%Acoeff(:,:,:,155) = c(5,6,:,:,:) *zetax

           Tdomain%specel(n)%Acoeff(:,:,:,156) = c(6,6,:,:,:) *xiz
           Tdomain%specel(n)%Acoeff(:,:,:,157) = c(6,6,:,:,:) *etaz
           Tdomain%specel(n)%Acoeff(:,:,:,158) = c(6,6,:,:,:) *zetaz

           Tdomain%specel(n)%Acoeff(:,:,:,159) = c(6,6,:,:,:) *xiy
           Tdomain%specel(n)%Acoeff(:,:,:,160) = c(6,6,:,:,:) *etay
           Tdomain%specel(n)%Acoeff(:,:,:,161) = c(6,6,:,:,:) *zetay
           !------------------------------------------------- 6

           !*********************************************************
           Tdomain%specel(n)%Acoeff(:,:,:,162) = -Whei * xix * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,163) = -Whei * xiy * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,164) = -Whei * xiz * Jac

           Tdomain%specel(n)%Acoeff(:,:,:,165) = -Whei * etax * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,166) = -Whei * etay * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,167) = -Whei * etaz * Jac           

           Tdomain%specel(n)%Acoeff(:,:,:,168) = -Whei * zetax * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,169) = -Whei * zetay * Jac
           Tdomain%specel(n)%Acoeff(:,:,:,170) = -Whei * zetaz * Jac
        endif
     endif

     !PML configuration
     !! Orientation
     if (Tdomain%specel(n)%PML) then
        allocate (wx (0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate (wy (0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate (wz (0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate (Id (0:ngllx-1,0:nglly-1,0:ngllz-1))

        if (Tdomain%sSubDomain(mat)%Px) then
           idef = Tdomain%specel(n)%Iglobnum(0,0,0)
           dx = Tdomain%GlobCoord(0,idef)
           idef = Tdomain%specel(n)%Iglobnum(ngllx-1,0,0)
           dx = abs(Tdomain%GlobCoord(0,idef) - dx) 
           if (Tdomain%sSubDomain(mat)%Left) then
              do i = 0,ngllx-1
                 ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(ngllx-1-i)) * float(ngllx-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(1,1,i,0,0)/Tdomain%specel(n)%Density(i,0,0))
                 else
                    vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
                    vp = sqrt(vp)
                 endif
                 wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wx = 0.
              enddo
           else 
              do i = 0,ngllx-1
                 ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(i)) * float(ngllx-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(1,1,i,0,0)/Tdomain%specel(n)%Density(i,0,0))
                 else
                    vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
                    vp = sqrt(vp)
                 endif
                 wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wx = 0.
              enddo
           endif
        else 
           wx = 0.

        endif

        if (Tdomain%sSubDomain(mat)%Py) then
           idef = Tdomain%specel(n)%Iglobnum(0,0,0)
           dx = Tdomain%GlobCoord(1,idef)
           idef = Tdomain%specel(n)%Iglobnum(0,nglly-1,0)
           dx = abs(Tdomain%GlobCoord (1,idef) - dx) 
           if (Tdomain%sSubDomain(mat)%Forward) then
              do j = 0,nglly-1
                 rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(nglly-1-j)) * float(nglly-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(2,2,0,j,0)/Tdomain%specel(n)%Density(0,j,0))
                 else
                    vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
                    vp = sqrt(vp)

                 endif
                 wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wy = 0.
              enddo
           else 
              do j = 0,nglly-1
                 rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(j)) * float(nglly-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(2,2,0,j,0)/Tdomain%specel(n)%Density(0,j,0))
                 else
                    vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
                    vp = sqrt(vp)
                 endif
                 wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wy = 0.
              enddo
           endif
        else 
           wy = 0.

        endif

        if (Tdomain%sSubDomain(mat)%Pz) then
           idef = Tdomain%specel(n)%Iglobnum(0,0,0)
           dx = Tdomain%GlobCoord(2,idef)
           idef = Tdomain%specel(n)%Iglobnum(0,0,ngllz-1)
           dx = abs(Tdomain%GlobCoord(2,idef) - dx)
           if (Tdomain%sSubDomain(mat)%Down) then
              do k = 0,ngllz-1
                 rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(ngllz-1-k)) * float(ngllz-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(3,3,0,0,k)/Tdomain%specel(n)%Density(0,0,k))
                 else
                    vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
                    vp = sqrt(vp)
                 endif
                 wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wz = 0.
              enddo
           else 
              do k = 0,ngllz-1
                 rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(k)) * float(ngllz-1)
                 if (Tdomain%specel(n)%aniso) then
                    vp=sqrt(c(3,3,0,0,k)/Tdomain%specel(n)%Density(0,0,k))
                 else
                    vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
                    vp = sqrt(vp)

                 endif
                 wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                      Tdomain%sSubdomain(mat)%npow)
                 !wz = 0.
              enddo
           endif
        else 
           wz = 0.

        endif
        
!!$        print *, 'define_array::coPMLrate', Tdomain%sSubdomain(mat)%coPMLrate       
        do i=0,ngllx-1
           do j=0,nglly-1
              do k=0,ngllz-1
                 maxRate=max(wx(i,j,k),wy(i,j,k),wz(i,j,k))
!!$                 print *, 'define_array::maxrate',wx(i,j,k),wy(i,j,k),wz(i,j,k),maxRate

                 if (wx(i,j,k)==0) wx(i,j,k)=Tdomain%sSubdomain(mat)%coPMLrate*maxRate
                 if (wy(i,j,k)==0) wy(i,j,k)=Tdomain%sSubdomain(mat)%coPMLrate*maxRate
                 if (wz(i,j,k)==0) wz(i,j,k)=Tdomain%sSubdomain(mat)%coPMLrate*maxRate

              enddo
           enddo
        enddo

        Id = 1.

        !! Strong formulation for stresses with respect to PML or Filtering PML
        if (Tdomain%specel(n)%FPML) then !Filtering PML

           Tdomain%specel(n)%DumpSx(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wx * Tdomain%sSubdomain(mat)%freq 
           Tdomain%specel(n)%DumpSx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,:,1)
           Tdomain%specel(n)%DumpSx (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wx * Tdomain%sSubdomain(mat)%freq) * &
                Tdomain%specel(n)%DumpSx(:,:,:,1)

           Tdomain%specel(n)%DumpSy(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wy * Tdomain%sSubdomain(mat)%freq 
           Tdomain%specel(n)%DumpSy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSy (:,:,:,1)
           Tdomain%specel(n)%DumpSy (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wy * Tdomain%sSubdomain(mat)%freq) * &
                Tdomain%specel(n)%DumpSy(:,:,:,1)

           Tdomain%specel(n)%DumpSz(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wz * Tdomain%sSubdomain(mat)%freq 
           Tdomain%specel(n)%DumpSz (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSz (:,:,:,1)
           Tdomain%specel(n)%DumpSz (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wz * Tdomain%sSubdomain(mat)%freq) *  &
                Tdomain%specel(n)%DumpSz(:,:,:,1)

           Tdomain%specel(n)%DumpMass(:,:,:,0) =  Tdomain%specel(n)%Density * Whei * Jac * wx * & 
                Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq
           Tdomain%specel(n)%DumpMass(:,:,:,1) =  Tdomain%specel(n)%Density * Whei * Jac * wy * & 
                Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq
           Tdomain%specel(n)%DumpMass(:,:,:,2) =  Tdomain%specel(n)%Density * Whei * Jac * wz * & 
                Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq

           Tdomain%specel(n)%Isx = Tdomain%sSubdomain(mat)%Dt * wx * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSx(:,:,:,1)     
           Tdomain%specel(n)%Isy = Tdomain%sSubdomain(mat)%Dt * wy * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSy(:,:,:,1)     
           Tdomain%specel(n)%Isz = Tdomain%sSubdomain(mat)%Dt * wz * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSz(:,:,:,1)     

           Tdomain%specel(n)%Ivx = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wx * Jac * Tdomain%sSubdomain(mat)%freq 
           Tdomain%specel(n)%Ivy = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wy * Jac * Tdomain%sSubdomain(mat)%freq 
           Tdomain%specel(n)%Ivz = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wz * Jac * Tdomain%sSubdomain(mat)%freq 

        else ! Normal PML

           Tdomain%specel(n)%DumpSx(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wx
           Tdomain%specel(n)%DumpSx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,:,1)
           Tdomain%specel(n)%DumpSx (:,:,:,0) = (Id- Tdomain%sSubdomain(mat)%Dt * 0.5 * wx) * Tdomain%specel(n)%DumpSx(:,:,:,1)     

           Tdomain%specel(n)%DumpSy(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wy
           Tdomain%specel(n)%DumpSy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSy (:,:,:,1) 
           Tdomain%specel(n)%DumpSy (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wy) * Tdomain%specel(n)%DumpSy(:,:,:,1)

           Tdomain%specel(n)%DumpSz(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wz
           Tdomain%specel(n)%DumpSz(:,:,:,1)  = 1./ Tdomain%specel(n)%DumpSz (:,:,:,1) 
           Tdomain%specel(n)%DumpSz (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wz) * Tdomain%specel(n)%DumpSz(:,:,:,1)

           Tdomain%specel(n)%DumpMass(:,:,:,0) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                Tdomain%sSubdomain(mat)%Dt * wx * Jac
           Tdomain%specel(n)%DumpMass(:,:,:,1) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                Tdomain%sSubdomain(mat)%Dt * wy * Jac
           Tdomain%specel(n)%DumpMass(:,:,:,2) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                Tdomain%sSubdomain(mat)%Dt * wz * Jac
        endif
        deallocate (wx,wy,wz,Id)
     endif

     deallocate (Jac, xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz, Whei)
     if (Tdomain%specel(n)%aniso) then
        deallocate(c)
     else
        deallocate(RKmod, Rmu, Rlam)
     endif
  enddo
!print *,'define_arrayyyyyyyyyy'
  !---------------------------------------------------

  ! Mass and DumpMass Communications inside Processors
  do n = 0,Tdomain%n_elem-1
     call get_Mass_Elem2Face(Tdomain,n)
     call get_Mass_Elem2Edge(Tdomain,n)
     call get_Mass_Elem2Vertex(Tdomain,n)
  enddo

  ! Invert Mass Matrix expression
  do n = 0,Tdomain%n_elem-1

     ngllx = Tdomain%specel(n)%ngllx
     nglly = Tdomain%specel(n)%nglly
     ngllz = Tdomain%specel(n)%ngllz
     allocate (LocMassMat(1:ngllx-2,1:nglly-2,1:ngllz-2))
     LocMassMat(:,:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)

     if (Tdomain%specel(n)%PML) then

        Tdomain%specel(n)%DumpVx (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
        Tdomain%specel(n)%DumpVx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVx (:,:,:,1) 
        Tdomain%specel(n)%DumpVx (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
        Tdomain%specel(n)%DumpVx (:,:,:,0) = Tdomain%specel(n)%DumpVx (:,:,:,0) * Tdomain%specel(n)%DumpVx (:,:,:,1)

        Tdomain%specel(n)%DumpVy (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
        Tdomain%specel(n)%DumpVy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVy (:,:,:,1) 
        Tdomain%specel(n)%DumpVy (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
        Tdomain%specel(n)%DumpVy (:,:,:,0) = Tdomain%specel(n)%DumpVy (:,:,:,0) * Tdomain%specel(n)%DumpVy (:,:,:,1)

        Tdomain%specel(n)%DumpVz (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
        Tdomain%specel(n)%DumpVz (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVz (:,:,:,1) 
        Tdomain%specel(n)%DumpVz (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
        Tdomain%specel(n)%DumpVz (:,:,:,0) = Tdomain%specel(n)%DumpVz (:,:,:,0) * Tdomain%specel(n)%DumpVz (:,:,:,1)

        if (Tdomain%specel(n)%FPML) then
           LocMassMat = Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,1:ngllz-2)
           deallocate (Tdomain%specel(n)%Ivx)
           allocate (Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,1:ngllz-2) )
           Tdomain%specel(n)%Ivx = LocMassMat *  Tdomain%specel(n)%DumpVx (:,:,:,1)
           LocMassMat = Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,1:ngllz-2)
           deallocate (Tdomain%specel(n)%Ivy)
           allocate (Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,1:ngllz-2) )
           Tdomain%specel(n)%Ivy = LocMassMat *  Tdomain%specel(n)%DumpVy (:,:,:,1)
           LocMassMat = Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,1:ngllz-2)
           deallocate (Tdomain%specel(n)%Ivz)
           allocate (Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,1:ngllz-2) )
           Tdomain%specel(n)%Ivz = LocMassMat *  Tdomain%specel(n)%DumpVz (:,:,:,1)
        endif

        deallocate (Tdomain%specel(n)%DumpMass)

     endif

     LocMassmat = 1./ LocMassMat
     deallocate (Tdomain%specel(n)%MassMat) 
     allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2) )
     Tdomain%specel(n)%MassMat = LocMassMat
     deallocate (LocMassMat)

     if (.not. Tdomain%specel(n)%aniso) then
        deallocate (Tdomain%specel(n)%Lambda)
        deallocate (Tdomain%specel(n)%Mu)  
     endif
!!$     if ((.not. Tdomain%LogicD%save_snapshots) .and. (.not. Tdomain%LogicD%save_energy)) then
!!$        deallocate (Tdomain%specel(n)%InvGrad)
!!$     else
!!$        if (((.not. Tdomain%Field_Order(2)).and. (.not. Tdomain%Field_Order(7))) .or. (Tdomain%specel(n)%PML) ) then
!!$           deallocate (Tdomain%specel(n)%InvGrad)
!!$        endif
!!$     endif
  enddo

  ! Define super objects properties (Btn)
  if (Tdomain%logicD%super_object_local_present) then
     if (Tdomain%super_object_type == "P") then

        do nf = 0, Tdomain%sPlaneW%n_faces-1
           ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
           ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2
           mat = Tdomain%sPlaneW%pFace(nf)%mat_index

           if (Tdomain%sPlaneW%pFace(nf)%dir == 0 .or. Tdomain%sPlaneW%pFace(nf)%dir == 5 ) then
              do j = 0,ngll2-1
                 do i = 0,ngll1-1
                    Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j)* &
                         Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
                 enddo
              enddo
           else if (Tdomain%sPlaneW%pFace(nf)%dir ==1 .or. Tdomain%sPlaneW%pFace(nf)%dir ==3 ) then
              do j = 0,ngll2-1
                 do i = 0,ngll1-1
                    Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                         Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
                 enddo
              enddo
           else
              do j = 0,ngll2-1
                 do i = 0,ngll1-1
                    Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwy(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                         Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
                 enddo
              enddo
           endif

           ! Internal communication of Btn
           do i = 0,3   
              ne = Tdomain%sPlaneW%pFace(nf)%Near_Edges(i)
              ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
              if ( Tdomain%sPlaneW%pFace(nf)%Orient_Edges(i) == 0 ) then
                 select case (i)
                 case (0)          
                    Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                         Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,0,0:2)
                 case (1)          
                    Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                         Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,1:ngll2-2,0:2)
                 case (2)          
                    Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                         Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,ngll2-1,0:2)
                 case (3)          
                    Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                         Tdomain%sPlaneW%pFace(nf)%Btn(0,1:ngll2-2,0:2)
                 end select
              else  
                 select case (i)
                 case (0)         
                    do j=1,ngll-2
                       Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1-j,0,0:2)
                    enddo
                 case (1) 
                    do j=1,ngll-2
                       Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,ngll2-1-j,0:2)
                    enddo
                 case (2)
                    do j=1,ngll-2          
                       Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1-j,ngll2-1,0:2)
                    enddo
                 case (3)
                    do j=1,ngll-2          
                       Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                            Tdomain%sPlaneW%pFace(nf)%Btn(0,ngll2-1-j,0:2)
                    enddo
                 end select
              endif

              nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(i)
              select case (i)
              case (0)          
                 Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                      Tdomain%sPlaneW%pFace(nf)%Btn(0,0,0:2)
              case (1)          
                 Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                      Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,0,0:2)
              case (2)          
                 Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                      Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,ngll2-1,0:2)
              case (3)          
                 Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                      Tdomain%sPlaneW%pFace(nf)%Btn(0,ngll2-1,0:2)
              end select
           enddo

           allocate (Store_Btn(1:ngll1-2,1:ngll2-2,0:2))
           Store_Btn (1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2)
           deallocate(Tdomain%sPlaneW%pFace(nf)%Btn)
           allocate(Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2))
           Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2) = Store_Btn(1:ngll1-2,1:ngll2-2,0:2)
           deallocate (Store_Btn)
           deallocate (Tdomain%sPlaneW%pFace(nf)%normal)

        enddo

     endif
  endif


  ! Define Neumann properties (Btn)
  if (Tdomain%logicD%neumann_local_present) then

     do nf = 0, Tdomain%sNeu%n_faces-1

        ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
        ngll2 = Tdomain%sNeu%nFace(nf)%ngll2
        mat = Tdomain%sNeu%nFace(nf)%mat_index

        if (Tdomain%sNeu%nFace(nf)%dir == 0 .or. Tdomain%sNeu%nFace(nf)%dir == 5 ) then
           do j = 0,ngll2-1
              do i = 0,ngll1-1
                 Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j)* &
                      Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
              enddo
           enddo
        else if (Tdomain%sNeu%nFace(nf)%dir ==1 .or. Tdomain%sNeu%nFace(nf)%dir ==3 ) then
           do j = 0,ngll2-1
              do i = 0,ngll1-1
                 Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                      Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
              enddo
           enddo
        else
           do j = 0,ngll2-1
              do i = 0,ngll1-1
                 Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwy(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                      Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
              enddo
           enddo
        endif

        ! Internal communication of Btn
        do i = 0,3   
           ne = Tdomain%sNeu%nFace(nf)%Near_Edges(i)
           ngll = Tdomain%sNeu%nEdge(ne)%ngll
           if ( Tdomain%sNeu%nFace(nf)%Orient_Edges(i) == 0 ) then
              select case (i)
              case (0)          
                 Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                      Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,0,0:2)
              case (1)          
                 Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                      Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,1:ngll2-2,0:2)
              case (2)          
                 Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                      Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,ngll2-1,0:2)
              case (3)          
                 Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                      Tdomain%sNeu%nFace(nf)%Btn(0,1:ngll2-2,0:2)
              end select
           else  
              select case (i)
              case (0)         
                 do j=1,ngll-2
                    Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                         Tdomain%sNeu%nFace(nf)%Btn(ngll1-1-j,0,0:2)
                 enddo
              case (1) 
                 do j=1,ngll-2
                    Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                         Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,ngll2-1-j,0:2)
                 enddo
              case (2)
                 do j=1,ngll-2          
                    Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                         Tdomain%sNeu%nFace(nf)%Btn(ngll1-1-j,ngll2-1,0:2)
                 enddo
              case (3)
                 do j=1,ngll-2          
                    Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                         Tdomain%sNeu%nFace(nf)%Btn(0,ngll2-1-j,0:2)
                 enddo
              end select
           endif

           nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(i)
           select case (i)
           case (0)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                   Tdomain%sNeu%nFace(nf)%Btn(0,0,0:2)
           case (1)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                   Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,0,0:2)
           case (2)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                   Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,ngll2-1,0:2)
           case (3)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                   Tdomain%sNeu%nFace(nf)%Btn(0,ngll2-1,0:2)
           end select
        enddo

        allocate (Store_Btn(1:ngll1-2,1:ngll2-2,0:2))
        Store_Btn (1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2)
        deallocate(Tdomain%sNeu%nFace(nf)%Btn)
        allocate(Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2))
        Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2) = Store_Btn(1:ngll1-2,1:ngll2-2,0:2)
        deallocate (Store_Btn)
        deallocate (Tdomain%sNeu%nFace(nf)%normal)

     enddo
  endif


  ! MPI communications
  if ( Tdomain%n_proc > 1 ) then

     do n = 0,Tdomain%n_proc-1

        ngll = 0
        ngllPML = 0
        ngllSO = 0

        do i = 0,Tdomain%sComm(n)%nb_faces-1

           nf = Tdomain%sComm(n)%faces(i)
           do j = 1,Tdomain%sFace(nf)%ngll2-2
              do k = 1,Tdomain%sFace(nf)%ngll1-2
                 Tdomain%sComm(n)%Give(ngll) = Tdomain%sFace(nf)%MassMat(k,j)
                 ngll = ngll + 1

              enddo
           enddo

           if (Tdomain%sFace(nf)%PML) then
              do j = 1,Tdomain%sFace(nf)%ngll2-2
                 do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2)
                    if (Tdomain%any_FPML) then
                       Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sFace(nf)%Ivx(k,j)
                       Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sFace(nf)%Ivy(k,j)
                       Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sFace(nf)%Ivz(k,j)

                    endif
                    ngllPML = ngllPML + 1

                 enddo
              enddo
           endif
        enddo

        do i = 0,Tdomain%sComm(n)%nb_edges-1
           ne = Tdomain%sComm(n)%edges(i)
           do j = 1,Tdomain%sEdge(ne)%ngll-2
              Tdomain%sComm(n)%Give(ngll) = Tdomain%sEdge(ne)%MassMat(j)
              ngll = ngll + 1
           enddo
           if (Tdomain%sEdge(ne)%PML) then
              do j = 1,Tdomain%sEdge(ne)%ngll-2
                 Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2)
                 if (Tdomain%any_FPML) then
                    Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sEdge(ne)%Ivx(j)
                    Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sEdge(ne)%Ivy(j)
                    Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sEdge(ne)%Ivz(j)
                 endif
                 ngllPML = ngllPML + 1
              enddo
           endif
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices-1
           nv =  Tdomain%sComm(n)%vertices(i)
           Tdomain%sComm(n)%Give(ngll) = Tdomain%svertex(nv)%MassMat
           ngll = ngll + 1
           if (Tdomain%sVertex(nv)%PML) then
              Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sVertex(nv)%DumpMass(0:2)
              if (Tdomain%any_FPML) then
                 Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sVertex(nv)%Ivx(0)
                 Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sVertex(nv)%Ivy(0)
                 Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sVertex(nv)%Ivz(0)
              endif
              ngllPML = ngllPML + 1
           endif
        enddo
        do i = 0,Tdomain%sComm(n)%nb_edges_so-1
           ne = Tdomain%sComm(n)%edges_SO(i)
           do j = 1,Tdomain%sPlaneW%pEdge(ne)%ngll-2
              Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2)
              ngllSO = ngllSO + 1
           enddo
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
           nv = Tdomain%sComm(n)%vertices_SO(i)
           Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2)
           ngllSO = ngllSO + 1    
        enddo
        do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
           ne = Tdomain%sComm(n)%edges_Neu(i)
           do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
              Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2)
              ngllSO = ngllSO + 1
           enddo
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
           nv = Tdomain%sComm(n)%vertices_Neu(i)
           Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2)
           ngllSO = ngllSO + 1    
        enddo
     enddo


     !do n = 0,Tdomain%n_proc-1
     !    if (Tdomain%sComm(n)%ngll>0) then
     !        call MPI_SEND (Tdomain%sComm(n)%Give, Tdomain%sComm(n)%ngll, MPI_DOUBLE_PRECISION, &
     !                       n, etiquette, MPI_COMM_WORLD, code)
     !        call MPI_RECV (Tdomain%sComm(n)%Take, Tdomain%sComm(n)%ngll, MPI_DOUBLE_PRECISION, &
     !                       n, etiquette, MPI_COMM_WORLD, statut, code)
     !    endif
     !    if (Tdomain%sComm(n)%ngllPML>0) then
     !        call MPI_SEND (Tdomain%sComm(n)%GivePML, 3*Tdomain%sComm(n)%ngllPML, MPI_DOUBLE_PRECISION, &
     !                       n, etiquette, MPI_COMM_WORLD, code)
     !        call MPI_RECV (Tdomain%sComm(n)%TakePML, 3*Tdomain%sComm(n)%ngllPML, MPI_DOUBLE_PRECISION, &
     !                       n, etiquette, MPI_COMM_WORLD, statut, code)
     !    endif
     !enddo


     n = Tdomain%n_proc
     do shift = 1,n-1
        I_give_to = rg + shift
        if (I_give_to > n-1)   I_give_to = I_give_to - n
        I_take_from = rg - shift
        if (I_take_from < 0)   I_take_from = I_take_from + n
        if (mod(n,shift)==0 .and. shift/=1) then
           n_rings = shift
        else if (mod(n,n-shift)==0 .and. shift/=n-1) then
           n_rings = n-shift
        else if (mod(n,2)==0 .and. mod(shift,2)==0) then
           n_rings = 2
        else
           n_rings = 1
        endif
        do i = 0,n_rings-1
           if (rg==i) then
              if (Tdomain%sComm(I_give_to)%ngll>0) then
                 call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                      MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
              endif
              if (Tdomain%sComm(I_take_from)%ngll>0) then
                 call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                      MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
              endif
           else
              do j = 0,n/n_rings-1
                 if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngll>0) then
                       call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngll>0) then
                       call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                 endif
              enddo
           endif
        enddo
        ! Peut-etre ce MPI_BARRIER n'est-il pas necessaire ?
        call MPI_BARRIER (MPI_COMM_WORLD, code)
        do i = 0,n_rings-1
           if (rg==i) then
              if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                 if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                    if (Tdomain%any_FPML) then
                       call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 6*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    else
                       call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                 endif
              endif
              if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                 if (Tdomain%any_FPML) then
                    call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 6*Tdomain%sComm(I_take_from)%ngllPML, &
                         MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                 else
                    call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                         MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                 endif
              endif
           else
              do j = 0,n/n_rings-1
                 if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                       if (Tdomain%any_FPML) then
                          call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 6*Tdomain%sComm(I_take_from)%ngllPML, &
                               MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                       else
                          call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                               MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                       endif
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                       if (Tdomain%any_FPML) then
                          call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 6*Tdomain%sComm(I_give_to)%ngllPML, &
                               MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                       else
                          call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                               MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                       endif
                    endif
                 endif
              enddo
           endif
        enddo
        ! Peut-etre ce MPI_BARRIER n'est-il pas necessaire ?
        call MPI_BARRIER (MPI_COMM_WORLD, code)
        do i = 0,n_rings-1
           if (rg==i) then
              if (Tdomain%sComm(I_give_to)%ngllSO>0) then
                 call MPI_SEND (Tdomain%sComm(I_give_to)%GiveSO, 3*Tdomain%sComm(I_give_to)%ngllSO, &
                      MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
              endif
              if (Tdomain%sComm(I_take_from)%ngllSO>0) then
                 call MPI_RECV (Tdomain%sComm(I_take_from)%TakeSO, 3*Tdomain%sComm(I_take_from)%ngllSO, &
                      MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
              endif
           else
              do j = 0,n/n_rings-1
                 if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllSO>0) then
                       call MPI_RECV (Tdomain%sComm(I_take_from)%TakeSO, 3*Tdomain%sComm(I_take_from)%ngllSO, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllSO>0) then
                       call MPI_SEND (Tdomain%sComm(I_give_to)%GiveSO, 3*Tdomain%sComm(I_give_to)%ngllSO, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                 endif
              enddo
           endif
        enddo
        ! Peut-etre ce MPI_BARRIER n'est-il pas necessaire ?
        call MPI_BARRIER (MPI_COMM_WORLD, code)
     enddo

     do n = 0,Tdomain%n_proc-1
        ngll = 0
        ngllPML = 0
        ngllSO = 0
        call Comm_Mass_Face(Tdomain,n,ngll,ngllPML)
        call Comm_Mass_Edge(Tdomain,n,ngll,ngllPML)
        call Comm_Mass_Vertex(Tdomain,n,ngll,ngllPML,rg)

        ! Super Object  (PlaneW or Fault) 
        do i = 0,Tdomain%sComm(n)%nb_edges_so-1
           ne = Tdomain%sComm(n)%edges_SO(i) 
           ngll1 = Tdomain%sPlaneW%pEdge(ne)%ngll
           if ( Tdomain%sComm(n)%orient_edges_SO(i) == 0 ) then 
              do j = 1,Tdomain%sPlaneW%pEdge(ne)%ngll-2
                 Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%splaneW%pEdge(ne)%Btn(j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
                 ngllSO = ngllSO + 1
              enddo
           else if ( Tdomain%sComm(n)%orient_edges_SO(i) == 1 ) then 
              do j = 1,Tdomain%splaneW%pEdge(ne)%ngll-2
                 Tdomain%sPlaneW%pEdge(ne)%Btn(ngll1-1-j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
                 ngllSO = ngllSO + 1
              enddo
           else 
              print*,'Pb with coherency number for edge' 
           endif
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
           nv = Tdomain%sComm(n)%vertices_SO(i) 
           Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%splaneW%pVertex(nv)%Btn(0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
           ngllSO = ngllSO + 1
        enddo

        ! Neumann
        do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
           ne = Tdomain%sComm(n)%edges_Neu(i) 
           ngll1 = Tdomain%sNeu%nEdge(ne)%ngll
           if ( Tdomain%sComm(n)%orient_edges_Neu(i) == 0 ) then 
              do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
                 Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
                 ngllSO = ngllSO + 1
              enddo
           else if ( Tdomain%sComm(n)%orient_edges_Neu(i) == 1 ) then 
              do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
                 Tdomain%sNeu%nEdge(ne)%Btn(ngll1-1-j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
                 ngllSO = ngllSO + 1
              enddo
           else 
              print*,'Pb with coherency number for edge' 
           endif
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
           nv = Tdomain%sComm(n)%vertices_Neu(i) 
           Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
           ngllSO = ngllSO + 1
        enddo

!!! A make "OPT = -O0 -C" compilation doesn't like these following deallocate. I don't know why !?!
        if (Tdomain%sComm(n)%ngll>0) then
           deallocate (Tdomain%sComm(n)%Give)
           deallocate (Tdomain%sComm(n)%Take)
        endif
        if (Tdomain%sComm(n)%ngllPML>0) then
           deallocate (Tdomain%sComm(n)%GivePML)
           deallocate (Tdomain%sComm(n)%TakePML)
        endif
        if (Tdomain%sComm(n)%ngllSO>0) then
           deallocate (Tdomain%sComm(n)%GiveSO)
           deallocate (Tdomain%sComm(n)%TakeSO)
        endif

     enddo
  endif


  do nf = 0,Tdomain%n_face-1
     if (Tdomain%sFace(nf)%PML) then
        Tdomain%sFace(nf)%DumpVx(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,1) = 1./Tdomain%sFace(nf)%DumpVx(:,:,1)
        Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%DumpVx(:,:,0) * Tdomain%sFace(nf)%DumpVx(:,:,1)     

        Tdomain%sFace(nf)%DumpVy(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,1) = 1./Tdomain%sFace(nf)%DumpVy (:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%DumpVy(:,:,0) * Tdomain%sFace(nf)%DumpVy(:,:,1)         

        Tdomain%sFace(nf)%DumpVz(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,1) = 1./Tdomain%sFace(nf)%DumpVz(:,:,1)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%DumpVz(:,:,0) * Tdomain%sFace(nf)%DumpVz(:,:,1)  
        if (Tdomain%sFace(nf)%FPML) then
           Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) *  Tdomain%sFace(nf)%DumpVx(:,:,1)
           Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) *  Tdomain%sFace(nf)%DumpVy(:,:,1)
           Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) *  Tdomain%sFace(nf)%DumpVz(:,:,1)
        endif
        deallocate (Tdomain%sFace(nf)%DumpMass)
     endif
     Tdomain%sFace(nf)%MassMat = 1./ Tdomain%sFace(nf)%MassMat
  enddo

  do ne = 0,Tdomain%n_edge-1
     if (Tdomain%sEdge(ne)%PML) then
        Tdomain%sEdge(ne)%DumpVx(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,1) = 1./Tdomain%sEdge(ne)%DumpVx(:,1)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%DumpVx(:,0) * Tdomain%sEdge(ne)%DumpVx(:,1)

        Tdomain%sEdge(ne)%DumpVy(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,1) = 1./Tdomain%sEdge(ne)%DumpVy(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%DumpVy(:,0) * Tdomain%sEdge(ne)%DumpVy(:,1)

        Tdomain%sEdge(ne)%DumpVz(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,1) = 1./Tdomain%sEdge(ne)%DumpVz(:,1)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%DumpVz(:,0) * Tdomain%sEdge(ne)%DumpVz(:,1)
        if (Tdomain%sEdge(ne)%FPML) then
           Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) *  Tdomain%sEdge(ne)%DumpVx(:,1)
           Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) *  Tdomain%sEdge(ne)%DumpVy(:,1)
           Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) *  Tdomain%sEdge(ne)%DumpVz(:,1)
        endif
        deallocate (Tdomain%sEdge(ne)%DumpMass)
     endif
     Tdomain%sEdge(ne)%MassMat = 1./ Tdomain%sEdge(ne)%MassMat
  enddo

  do nv = 0,Tdomain%n_vertex-1
     if (Tdomain%sVertex(nv)%PML) then
        Tdomain%sVertex(nv)%DumpVx(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(1) = 1./Tdomain%sVertex(nv)%DumpVx(1)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%DumpVx(0) * Tdomain%sVertex(nv)%DumpVx(1)

        Tdomain%sVertex(nv)%DumpVy(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(1) = 1./Tdomain%sVertex(nv)%DumpVy(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%DumpVy(0) * Tdomain%sVertex(nv)%DumpVy(1)

        Tdomain%sVertex(nv)%DumpVz(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(1) = 1./Tdomain%sVertex(nv)%DumpVz(1)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%DumpVz(0) * Tdomain%sVertex(nv)%DumpVz(1)
        if (Tdomain%sVertex(nv)%FPML) then
           Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) *  Tdomain%sVertex(nv)%DumpVx(1)
           Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) *  Tdomain%sVertex(nv)%DumpVy(1)
           Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) *  Tdomain%sVertex(nv)%DumpVz(1)
        endif
     endif
     Tdomain%sVertex(nv)%MassMat = 1./ Tdomain%sVertex(nv)%MassMat
  enddo


  ! Mass Mat fo SO
  if (Tdomain%logicD%super_object_local_present) then
     if (Tdomain%super_object_type == "P") then

        do nf = 0, Tdomain%sPlaneW%n_faces-1     
           nv_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
           Tdomain%sPlaneW%pFace(nf)%MassMat_Up = Tdomain%sFace(nv_aus)%MassMat

           nv_aus = Tdomain%sPlaneW%pFace(nf)%Face_DOWN
           Tdomain%sPlaneW%pFace(nf)%MassMat_Down =  Tdomain%sFace(nv_aus)%MassMat

        enddo
        do ne = 0, Tdomain%sPlaneW%n_edges-1  
           nv_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
           Tdomain%sPlaneW%pEdge(ne)%MassMat_Up = Tdomain%sEdge(nv_aus)%MassMat

           nv_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN
           Tdomain%sPlaneW%pEdge(ne)%MassMat_Down =  Tdomain%sEdge(nv_aus)%MassMat

        enddo
        do nv = 0, Tdomain%sPlaneW%n_vertices-1
           nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
           Tdomain%sPlaneW%pVertex(nv)%MassMat_Up = Tdomain%sVertex(nv_aus)%MassMat

           nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
           Tdomain%sPlaneW%pVertex(nv)%MassMat_Down =  Tdomain%sVertex(nv_aus)%MassMat

        enddo

     endif
  endif
!print *, 'define_array::check0a',rg
  if (notOpen==0) call MPI_FILE_CLOSE(desc,code)
!print *, 'define_array::check1a',rg
  call system ('/bin/echo -e ''\a''')
  return
end subroutine Define_Arrays


! $Id: define_arrays.f90,v 1.10 2007/08/30 15:39:20 taquanganh Exp $ !
