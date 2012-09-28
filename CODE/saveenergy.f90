subroutine saveenergy(Tdomain,rg,icount) 

  ! Modified Regis 04/2007

  use sdomain

  implicit none

  include 'mpif.h'

  type (domain),intent (IN), target :: Tdomain
  integer, intent(IN) :: rg, icount
  real  :: Kin, Pot, Energy,KinGLL,PotGLL,zp,&
       Ecurl, Ediv, EcurlGLL, EdivGLL,PotP,PotS,PotResi !,PotP_GLL,PotS_GLL,PotResi_GLL

double precision :: PotP_GLL,PotS_GLL,PotResi_GLL


  integer :: n, i, j, k, n_z, n1, n2, n3,geo
  integer :: ngllx, nglly, ngllz, mat, desc, code, nbprocs, NbOctReal, NbOctInt
  integer (kind=MPI_OFFSET_KIND) :: PosFile
  integer, dimension(MPI_STATUS_SIZE) :: statut
  real, dimension(:,:,:), allocatable :: Whei, dUx_dxi, dUx_deta, dUx_dzeta, &
       dUy_dxi, dUy_deta, dUy_dzeta, &
       dUz_dxi, dUz_deta, dUz_dzeta,&
       dUx_dx, dUx_dy,dUx_dz,&
       dUy_dx, dUy_dy,dUy_dz,&
       dUz_dx, dUz_dy,dUz_dz,&
       xix ,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz
  real, dimension(:,:,:,:), allocatable :: Ediv_curl
  logical :: EnerCount

real :: kappa_eqv, mu_eqv
real, dimension(6,6) :: s_C

  Kin = 0.
  Pot = 0.
  Ecurl=0.
  Ediv=0.
  PotP=0.
  PotS=0.
  PotResi=0.

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,code)
  call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,code)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,"SaveEnergy.out",&
       MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)
  if ((icount==1).and.(rg==0)) then
     PosFile=0
     call MPI_FILE_WRITE_AT(desc,PosFile,nbprocs,1,MPI_INTEGER,statut,code)  
  endif
   !print *, 'SaveEner_check0'
  do n = 0, Tdomain%n_elem-1
  !print *, 'SaveEner_check1'
     mat = Tdomain%specel(n)%mat_index
     EnerCount=((mat==Tdomain%EnerMat-1) .or. (Tdomain%EnerMat==0))
     if (EnerCount) then 
     !if (mat==Tdomain%EnerMat-1) then 
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (.not. allocated (Whei)) allocate (Whei (0:ngllx-1,0:nglly-1,0:ngllz-1))
        do k = 0, ngllz -1 
           do j = 0, nglly -1 
              do i = 0, ngllx -1
                 Whei (i,j,k) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j) &
                      * Tdomain%sSubdomain(mat)%GLLwz(k)
              enddo
           enddo
        enddo
   !print *, 'SaveEner_check2'
        if (.not. Tdomain%specel(n)%PML) then

           call Displ_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
           call Displ_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
           call Displ_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)

!print *, 'SaveEner_check3'
           if (.not. allocated (dUx_dxi)) allocate (dUx_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_deta)) allocate (dUx_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_dzeta)) allocate (dUx_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dxi)) allocate (dUy_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_deta)) allocate (dUy_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dzeta)) allocate (dUy_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dxi)) allocate (dUz_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_deta)) allocate (dUz_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dzeta)) allocate (dUz_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
!print *, 'SaveEner_check4'
           call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                ngllx, Tdomain%specel(n)%Displ(0,0,0,0), ngllx, 0., dUx_dxi, ngllx)
           do n_z = 0,ngllz-1
              call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Displ(0,0,n_z,0),   &
                   ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dUx_deta(0,0,n_z), ngllx)
           enddo
           call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Displ(0,0,0,0),  &
                ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dUx_dzeta, ngllx*nglly)
           call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                ngllx,Tdomain%specel(n)%Displ(0,0,0,1), ngllx, 0., dUy_dxi, ngllx)
           do n_z = 0,ngllz-1
              call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Displ(0,0,n_z,1),   &
                   ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dUy_deta(0,0,n_z), ngllx)
           enddo
           call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Displ(0,0,0,1),   &
                ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dUy_dzeta, ngllx*nglly)
           call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,    &
                ngllx,Tdomain%specel(n)%Displ(0,0,0,2), ngllx, 0., dUz_dxi, ngllx)
           do n_z = 0,ngllz-1
              call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Displ(0,0,n_z,2),    &
                   ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dUz_deta(0,0,n_z), ngllx)
           enddo
           call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Displ(0,0,0,2),   &
                ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dUz_dzeta, ngllx*nglly)
!print *, 'SaveEner_check5'
!!$      n1=ngllx-2
!!$      n2=nglly-2
!!$      n3=ngllz-2
!!$      call DGEMM ( 'N', 'N', n1, n2*n3, n1, 1., Tdomain%sSubdomain(mat)%hTprimex(1,1),   &
!!$                        ngllx, Tdomain%specel(n)%Displ(1,1,1,0), n1, 0., dUx_dxi, n1)
!!$      do n_z = 1,n2
!!$         call DGEMM ( 'N', 'N', n1, n2, n2, 1., Tdomain%specel(n)%Displ(1,1,n_z,0),   &
!!$                           n1, Tdomain%sSubdomain(mat)%hprimey(1,1), nglly, 0., dUx_deta(1,1,n_z), n1)
!!$      enddo
!!$      call DGEMM ( 'N', 'N', n1*n2, n3, n3, 1., Tdomain%specel(n)%Displ(1,1,1,0),  &
!!$                        n1*n2, Tdomain%sSubdomain(mat)%hprimez(1,1), ngllz, 0., dUx_dzeta, n1*n2)
!!$      call DGEMM ( 'N', 'N', n1, n2*n3, n1, 1., Tdomain%sSubdomain(mat)%hTprimex(1,1),   &
!!$                        ngllx, Tdomain%specel(n)%Displ(1,1,1,1), n1, 0., dUy_dxi, n1)
!!$      do n_z = 1,n2
!!$         call DGEMM ( 'N', 'N', n1, n2, n2, 1., Tdomain%specel(n)%Displ(1,1,n_z,1),   &
!!$                           n1, Tdomain%sSubdomain(mat)%hprimey(1,1), nglly, 0., dUy_deta(1,1,n_z), n1)
!!$      enddo
!!$      call DGEMM ( 'N', 'N', n1*n2, n3, n3, 1., Tdomain%specel(n)%Displ(1,1,1,1),  &
!!$                        n1*n2, Tdomain%sSubdomain(mat)%hprimez(1,1), ngllz, 0., dUy_dzeta, n1*n2)
!!$      call DGEMM ( 'N', 'N', n1, n2*n3, n1, 1., Tdomain%sSubdomain(mat)%hTprimex(1,1),   &
!!$                        ngllx, Tdomain%specel(n)%Displ(1,1,1,2), n1, 0., dUz_dxi, n1)
!!$      do n_z = 1,n2
!!$         call DGEMM ( 'N', 'N', n1, n2, n2, 1., Tdomain%specel(n)%Displ(1,1,n_z,2),   &
!!$                           n1, Tdomain%sSubdomain(mat)%hprimey(1,1), nglly, 0., dUz_deta(1,1,n_z), n1)
!!$      enddo
!!$      call DGEMM ( 'N', 'N', n1*n2, n3, n3, 1., Tdomain%specel(n)%Displ(1,1,1,2),  &
!!$                        n1*n2, Tdomain%sSubdomain(mat)%hprimez(1,1), ngllz, 0., dUz_dzeta, n1*n2)


           call Veloc_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
           call Veloc_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
           call Veloc_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)
!print *, 'SaveEner_check6'

           !==========



           allocate (xix(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (xiy(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (xiz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etax(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etay(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (zetax(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (zetay(0:ngllx-1,0:nglly-1,0:ngllz-1)) 
           allocate (zetaz(0:ngllx-1,0:nglly-1,0:ngllz-1))

!print *, 'SaveEner_check7'
           xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)

           xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
           xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)

           etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
           etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
           etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)


           zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
           zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
           zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)

!print *, 'SaveEner_check8'
           if (.not. allocated (dUx_dx)) allocate (dUx_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_dy)) allocate (dUx_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_dz)) allocate (dUx_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dx)) allocate (dUy_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dy)) allocate (dUy_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dz)) allocate (dUy_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dx)) allocate (dUz_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dy)) allocate (dUz_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dz)) allocate (dUz_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))

!print *, 'SaveEner_check9'
           dUx_dx=dUx_dxi*xix+dUx_deta*etax+dUx_dzeta*zetax
           dUy_dx=dUy_dxi*xix+dUy_deta*etax+dUy_dzeta*zetax
           dUz_dx=dUz_dxi*xix+dUz_deta*etax+dUz_dzeta*zetax

           dUx_dy=dUx_dxi*xiy+dUx_deta*etay+dUx_dzeta*zetay
           dUy_dy=dUy_dxi*xiy+dUy_deta*etay+dUy_dzeta*zetay
           dUz_dy=dUz_dxi*xiy+dUz_deta*etay+dUz_dzeta*zetay

           dUx_dz=dUx_dxi*xiz+dUx_deta*etaz+dUx_dzeta*zetaz
           dUy_dz=dUy_dxi*xiz+dUy_deta*etaz+dUy_dzeta*zetaz
           dUz_dz=dUz_dxi*xiz+dUz_deta*etaz+dUz_dzeta*zetaz
           deallocate(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz)

           if (.not. allocated(Ediv_curl)) allocate(Ediv_curl(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
           Ediv_curl(:,:,:,0)=(dUx_dx+dUy_dy+dUz_dz)**2
           Ediv_curl(:,:,:,1)=(dUy_dz-dUz_dy)**2+(dUz_dx-dUx_dz)**2+(dUx_dy-dUy_dx)**2
           deallocate(dUx_dx,dUx_dy,dUx_dz,&
                dUy_dx,dUy_dy,dUy_dz,&
                dUz_dx,dUz_dy,dUz_dz)
           !==========



!print *, 'SaveEner_check9'







           do k = 0, ngllz-1 
              do j = 0, nglly-1 
                 do i = 0, ngllx-1



                    KinGLL=Tdomain%specel(n)%Density(i,j,k) * Whei(i,j,k) * Tdomain%specel(n)%Jacob(i,j,k) * &
                         (Tdomain%specel(n)%Veloc(i,j,k,0)**2+ Tdomain%specel(n)%Veloc(i,j,k,1)**2 + &
                         Tdomain%specel(n)%Veloc(i,j,k,2)**2)
                    PotGLL= -Tdomain%specel(n)%Acoeff(i,j,k,0) * (dUx_dxi(i,j,k))**2                 &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,1) * dUx_dxi(i,j,k)   * dUx_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,2) * dUx_dxi(i,j,k)   * dUx_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,3) * dUx_dxi(i,j,k)   * dUy_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,4) * dUx_dxi(i,j,k)   * dUy_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,5) * dUx_dxi(i,j,k)   * dUy_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,6) * dUx_dxi(i,j,k)   * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,7) * dUx_dxi(i,j,k)   * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,8) * dUx_dxi(i,j,k)   * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,9) * (dUx_deta(i,j,k))**2                &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,10)* dUx_deta(i,j,k)  * dUx_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,11)* dUx_deta(i,j,k)  * dUy_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,12)* dUx_deta(i,j,k)  * dUy_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,13)* dUx_deta(i,j,k)  * dUy_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,14)* dUx_deta(i,j,k)  * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,15)* dUx_deta(i,j,k)  * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,16)* dUx_deta(i,j,k)  * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,17)* (dUx_dzeta(i,j,k))**2               &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,18)* dUx_dzeta(i,j,k) * dUy_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,19)* dUx_dzeta(i,j,k) * dUy_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,20)* dUx_dzeta(i,j,k) * dUy_dzeta (i,j,k)&
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,21)* dUx_dzeta(i,j,k) * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,22)* dUx_dzeta(i,j,k) * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,23)* dUx_dzeta(i,j,k) * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,24)* (dUy_dxi(i,j,k))**2                 &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,25)* dUy_dxi(i,j,k)   * dUy_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,26)* dUy_dxi(i,j,k)   * dUy_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,27)* dUy_dxi(i,j,k)   * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,28)* dUy_dxi(i,j,k)   * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,29)* dUy_dxi(i,j,k)   * dUz_dzeta (i,j,k)&
                         -   Tdomain%specel(n)%Acoeff(i,j,k,30)* (dUy_deta(i,j,k))**2                &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,31)* dUy_deta(i,j,k)  * dUy_dzeta(i,j,k) &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,32)* dUy_deta(i,j,k)  * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,33)* dUy_deta(i,j,k)  * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,34)* dUy_deta(i,j,k)  * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,35)* (dUy_dzeta(i,j,k))**2               &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,36)* dUy_dzeta(i,j,k) * dUz_dxi(i,j,k)   &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,37)* dUy_dzeta(i,j,k) * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,38)* dUy_dzeta(i,j,k) * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,39)* (dUz_dxi(i,j,k))**2                 &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,40)* dUz_dxi(i,j,k)   * dUz_deta(i,j,k)  &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,41)* dUz_dxi(i,j,k)   * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,42)* (dUz_deta(i,j,k))**2                &
                         - 2*Tdomain%specel(n)%Acoeff(i,j,k,43)* dUz_deta(i,j,k)  * dUz_dzeta(i,j,k) &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,44)* (dUz_dzeta(i,j,k))**2

               s_C=Tdomain%specel(n)%c(:,:,i,j,k)
               kappa_eqv=(s_C(1,1)+s_C(2,2)+s_C(3,3)+2.*(s_C(1,2)+s_C(1,3)+s_C(2,3)))/9.
               mu_eqv=(s_C(1,1)+s_C(2,2)+s_C(3,3)+3*(s_C(4,4)+s_C(5,5)+s_C(6,6))-(s_C(1,2)+s_C(1,3)+s_C(2,3)))/15.
              !print *, 'kappa_eqv check', kappa_eqv

                    EdivGLL=Whei(i,j,k)*Ediv_curl(i,j,k,0)* Tdomain%specel(n)%Jacob(i,j,k)*(kappa_eqv+4*mu_eqv/3)/2
                    EcurlGLL=Whei(i,j,k)*Ediv_curl(i,j,k,1)* Tdomain%specel(n)%Jacob(i,j,k)*mu_eqv/2

                    !EdivGLL=Whei(i,j,k)*Ediv_curl(i,j,k,0)* Tdomain%specel(n)%Jacob(i,j,k)
                    !EcurlGLL=Whei(i,j,k)*Ediv_curl(i,j,k,1)* Tdomain%specel(n)%Jacob(i,j,k)


                    PotP_GLL=&
                         -Tdomain%specel(n)%Acoeff(i,j,k,0) * (dUx_dxi(i,j,k))**2   &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,30)* (dUy_deta(i,j,k))**2                &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,44)* (dUz_dzeta(i,j,k))**2 &
- 2*Tdomain%specel(n)%Acoeff(i,j,k,4) * dUx_dxi(i,j,k)   * dUy_deta(i,j,k)  &
 - 2*Tdomain%specel(n)%Acoeff(i,j,k,8) * dUx_dxi(i,j,k)   * dUz_dzeta(i,j,k) &
 - 2*Tdomain%specel(n)%Acoeff(i,j,k,34)* dUy_deta(i,j,k)  * dUz_dzeta(i,j,k) 
                    PotS_GLL=&
                         -   Tdomain%specel(n)%Acoeff(i,j,k,9) * (dUx_deta(i,j,k))**2                &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,24)* (dUy_dxi(i,j,k))**2                 &
                         + 2*Tdomain%specel(n)%Acoeff(i,j,k,11)* dUx_deta(i,j,k)  * dUy_dxi(i,j,k)   &
                         -  Tdomain%specel(n)%Acoeff(i,j,k,17)* (dUx_dzeta(i,j,k))**2               &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,39)* (dUz_dxi(i,j,k))**2                 &
                         +2*Tdomain%specel(n)%Acoeff(i,j,k,21)* dUx_dzeta(i,j,k) * dUz_dxi(i,j,k)   &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,35)* (dUy_dzeta(i,j,k))**2               &
                         -   Tdomain%specel(n)%Acoeff(i,j,k,42)* (dUz_deta(i,j,k))**2                &
                         +2*Tdomain%specel(n)%Acoeff(i,j,k,37)* dUy_dzeta(i,j,k) * dUz_deta(i,j,k)

                    !PotResi_GLL=PotGLL-PotP_GLL-PotS_GLL  
if (PotP_GLL<1e-200 .and. PotS_GLL<1e-200) then
PotP_GLL = 1e-200
PotS_GLL = 1e-200
endif

                    PotResi_GLL=PotP_GLL/PotS_GLL 

!print *, 'check_floating error', PotS_GLL, PotP_GLL

                    zp=Tdomain%GlobCoord (2,Tdomain%specel(n)%Iglobnum(i,j,k))    
                    geo=0
                    if (i==0 .or. i==ngllx-1) geo=geo+1
                    if (j==0 .or. j==nglly-1) geo=geo+1
                    if (k==0 .or. k==ngllz-1) geo=geo+1
                    if (zp==0) then
                       if (geo==0 .or. geo==1) then
                          Kin=Kin+KinGLL
                          Pot=Pot+PotGLL
                          Ediv=Ediv+EdivGLL
                          Ecurl=Ecurl+EcurlGLL 
                          PotP=PotP+PotP_GLL
                          PotS=PotS+PotS_GLL
                          PotResi=PotResi+PotResi_GLL
                       elseif (geo==2) then
                          Kin=Kin+KinGLL/2
                          Pot=Pot+PotGLL/2
                          Ediv=Ediv+EdivGLL/2
                          Ecurl=Ecurl+EcurlGLL/2
                          PotP=PotP+PotP_GLL/2
                          PotS=PotS+PotS_GLL/2
                          PotResi=PotResi+PotResi_GLL/2
                       elseif (geo==2) then
                          Kin=Kin+KinGLL/4
                          Pot=Pot+PotGLL/4
                          Ediv=Ediv+EdivGLL/4
                          Ecurl=Ecurl+EcurlGLL/4
                          PotP=PotP+PotP_GLL/4
                          PotS=PotS+PotS_GLL/4
                          PotResi=PotResi+PotResi_GLL/4
                       endif
                    else
                       if (geo==0) then
                          Kin=Kin+KinGLL
                          Pot=Pot+PotGLL
                          Ediv=Ediv+EdivGLL
                          Ecurl=Ecurl+EcurlGLL
                          PotP=PotP+PotP_GLL
                          PotS=PotS+PotS_GLL
                          PotResi=PotResi+PotResi_GLL
                       elseif (geo==1) then
                          Kin=Kin+KinGLL/2
                          Pot=Pot+PotGLL/2
                          Ediv=Ediv+EdivGLL/2
                          Ecurl=Ecurl+EcurlGLL/2
                          PotP=PotP+PotP_GLL/2
                          PotS=PotS+PotS_GLL/2
                          PotResi=PotResi+PotResi_GLL/2
                       elseif (geo==2) then
                          Kin=Kin+KinGLL/4
                          Pot=Pot+PotGLL/4
                          Ediv=Ediv+EdivGLL/4
                          Ecurl=Ecurl+EcurlGLL/4
                          PotP=PotP+PotP_GLL/4
                          PotS=PotS+PotS_GLL/4
                          PotResi=PotResi+PotResi_GLL/4
                       elseif (geo==2) then
                          Kin=Kin+KinGLL/8
                          Pot=Pot+PotGLL/8
                          Ediv=Ediv+EdivGLL/8
                          Ecurl=Ecurl+EcurlGLL/8
                          PotP=PotP+PotP_GLL/8
                          PotS=PotS+PotS_GLL/8
                          PotResi=PotResi+PotResi_GLL/8
                       endif
                    endif

                 enddo
              enddo
           enddo
        endif
     endif
  end do

  deallocate(Ediv_curl)
!!$Kin=0.5*Kin
!!$Pot=0.5*Pot
!!$Energy = Kin+Pot
!!$PosFile = NbOctInt + ( nbprocs*(icount-1) + rg )* NbOctReal*3
!!$call MPI_FILE_WRITE_AT(desc,PosFile,Kin,1,MPI_DOUBLE_PRECISION,statut,code)
!!$PosFile = PosFile + NbOctReal
!!$call MPI_FILE_WRITE_AT(desc,PosFile,Pot,1,MPI_DOUBLE_PRECISION,statut,code)
!!$PosFile = PosFile + NbOctReal
!!$call MPI_FILE_WRITE_AT(desc,PosFile,Energy,1,MPI_DOUBLE_PRECISION,statut,code)
!print *, 'SaveEner_checkA'

  Kin=0.5*Kin
  Pot=0.5*Pot
  PotP=0.5*PotP
  PotS=0.5*PotS
  PotResi=0.5*PotResi
  Energy = Kin+Pot
  PosFile = NbOctInt + ( nbprocs*(icount-1) + rg )* NbOctReal*8
  call MPI_FILE_WRITE_AT(desc,PosFile,Kin,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,Pot,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,Energy,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,Ediv,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,Ecurl,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,PotP,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,PotS,1,MPI_DOUBLE_PRECISION,statut,code)
  PosFile = PosFile + NbOctReal
  call MPI_FILE_WRITE_AT(desc,PosFile,PotResi,1,MPI_DOUBLE_PRECISION,statut,code)
!print *, 'SaveEner_checkB'

  call MPI_FILE_CLOSE(desc,code)

  return
end subroutine saveenergy





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from vertex to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displ_Collection_From_Vertex(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)


  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: nv, s_ngllx,s_nglly,s_ngllz
  nv = s_Tdomain%specel(s_n)%near_vertices(0)
  s_Tdomain%specel(s_n)%Displ(0,0,0,0) = s_Tdomain%sVertex(nv)%Displ(0) 
  s_Tdomain%specel(s_n)%Displ(0,0,0,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(0,0,0,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(1)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,0,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,0,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,0,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(2)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,0,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,0,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,0,2) = s_Tdomain%sVertex(nv)%Displ(2) 

  nv = s_Tdomain%specel(s_n)%near_vertices(3)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,0,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,0,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,0,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(4)
  s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(5)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(6)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Displ(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(7)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Displ(0)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Displ(1)
  s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Displ(2)
end subroutine Displ_Collection_From_Vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Face to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displ_Collection_From_Face(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)
  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: nf, s_ngllx,s_nglly,s_ngllz,ngll1,ngll2,i,j

  nf = s_Tdomain%specel(s_n)%near_faces(0)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,0,0) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,0,1) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,0,2) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,i,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,i,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,i,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,0,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,0,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,0,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(1)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(i,0,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(i,0,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(i,0,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,0,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,0,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,0,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,0,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,0,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,0,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)  
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,0,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(2)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Displ(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(3)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_ngllz-2,0) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_ngllz-2,1) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_ngllz-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(4)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Displ(1:s_nglly-2,1:s_ngllz-2,0) 
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Displ(1:s_nglly-2,1:s_ngllz-2,1) 
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Displ(1:s_nglly-2,1:s_ngllz-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,j,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,j,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,j,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0)
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1)
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2) 
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(5)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_nglly-2,0) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_nglly-2,1) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(1:s_ngllx-2,1:s_nglly-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(i,s_nglly-1-j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(j,s_nglly-1-i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Displ(i,j,0) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Displ(i,j,1) 
           s_Tdomain%specel(s_n)%Displ(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Displ(i,j,2)
        enddo
     enddo
  endif
end subroutine Displ_Collection_From_Face
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Edges to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displ_Collection_From_Edge(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)

  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: ne, s_ngllx,s_nglly,s_ngllz,ngll,i

  ne = s_Tdomain%specel(s_n)%near_edges(0)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(0) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,0,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0) 
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,0,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,0,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,0,0) = s_Tdomain%sEdge(ne)%Displ(i,0) 
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,0,1) = s_Tdomain%sEdge(ne)%Displ(i,1) 
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,0,2) = s_Tdomain%sEdge(ne)%Displ(i,2) 
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(1)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(1) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,0,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,0,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,0,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,0,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,0,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,0,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(2)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(2) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,0,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,0,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,0,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2) 
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,0,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,0,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,0,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(3)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(3) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,0,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,0,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,0,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,0,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,0,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,0,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(4)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(4) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,0,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(5)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(5) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,0,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,0,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(6)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(6) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(0,0,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(0,0,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(0,0,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(0,0,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(7)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(7) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(8)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(8) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(s_ngllx-1,1:s_nglly-2,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(i,0) 
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(i,1) 
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1,s_nglly-1-i,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(i,2) 
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(9)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(9) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(1:s_ngllx-2,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(s_ngllx-1-i,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(10)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(10) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(11)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(11) == 0 ) then
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Displ(0,1:s_nglly-2,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Displ(i,0)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Displ(i,1)
        s_Tdomain%specel(s_n)%Displ(0,s_nglly-1-i,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Displ(i,2)
     enddo
  endif
end subroutine Displ_Collection_From_Edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from vertex to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Veloc_Collection_From_Vertex(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)


  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: nv, s_ngllx,s_nglly,s_ngllz
  nv = s_Tdomain%specel(s_n)%near_vertices(0)
  s_Tdomain%specel(s_n)%Veloc(0,0,0,0) = s_Tdomain%sVertex(nv)%Veloc(0) 
  s_Tdomain%specel(s_n)%Veloc(0,0,0,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(0,0,0,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(1)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,0,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,0,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,0,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(2)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,0,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,0,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,0,2) = s_Tdomain%sVertex(nv)%Veloc(2) 

  nv = s_Tdomain%specel(s_n)%near_vertices(3)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,0,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,0,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,0,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(4)
  s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(5)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(6)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Veloc(2)

  nv = s_Tdomain%specel(s_n)%near_vertices(7)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sVertex(nv)%Veloc(0)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sVertex(nv)%Veloc(1)
  s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sVertex(nv)%Veloc(2)
end subroutine Veloc_Collection_From_Vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Face to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Veloc_Collection_From_Face(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)
  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: nf, s_ngllx,s_nglly,s_ngllz,ngll1,ngll2,i,j

  nf = s_Tdomain%specel(s_n)%near_faces(0)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,0,0) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,0,1) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,0,2) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,i,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,i,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,i,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(0) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,0,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,0,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,0,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(1)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(i,0,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(i,0,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(i,0,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,0,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,0,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,0,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,0,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,0,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,0,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)  
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(1) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,0,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(2)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,0) 
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,1) 
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Veloc(1:ngll1-2,1:ngll2-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(2) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(3)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_ngllz-2,0) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_ngllz-2,1) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_ngllz-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(3) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(4)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,1:s_ngllz-2,0) =  s_Tdomain%sFace(nf)%Veloc(1:s_nglly-2,1:s_ngllz-2,0) 
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,1:s_ngllz-2,1) =  s_Tdomain%sFace(nf)%Veloc(1:s_nglly-2,1:s_ngllz-2,1) 
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,1:s_ngllz-2,2) =  s_Tdomain%sFace(nf)%Veloc(1:s_nglly-2,1:s_ngllz-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1-j,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1-j,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1-j,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,j,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,j,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,j,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0)
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1)
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(4) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,s_ngllz-1-i,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,s_ngllz-1-i,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-j,s_ngllz-1-i,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2) 
        enddo
     enddo
  endif

  nf = s_Tdomain%specel(s_n)%near_faces(5)
  ngll1 = s_Tdomain%sFace(nf)%ngll1
  ngll2 = s_Tdomain%sFace(nf)%ngll2
  if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_nglly-2,0) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_nglly-2,1) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,1:s_nglly-2,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(1:s_ngllx-2,1:s_nglly-2,2) 
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 1 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 2 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(i,s_nglly-1-j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 3 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1-j,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 4 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 5 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 6 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(j,s_nglly-1-i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  else if ( s_Tdomain%specel(s_n)%Orient_Faces(5) == 7 ) then
     do i=1,ngll1-2
        do j=1,ngll2-2
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,0) =  s_Tdomain%sFace(nf)%Veloc(i,j,0) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,1) =  s_Tdomain%sFace(nf)%Veloc(i,j,1) 
           s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-j,s_nglly-1-i,s_ngllz-1,2) =  s_Tdomain%sFace(nf)%Veloc(i,j,2)
        enddo
     enddo
  endif
end subroutine Veloc_Collection_From_Face
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Edges to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Veloc_Collection_From_Edge(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)

  use sdomain

  implicit none

  type (Domain), intent (INOUT) :: s_Tdomain
  integer, intent (IN) :: s_n

  integer :: ne, s_ngllx,s_nglly,s_ngllz,ngll,i

  ne = s_Tdomain%specel(s_n)%near_edges(0)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(0) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,0,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0) 
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,0,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,0,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,0,0) = s_Tdomain%sEdge(ne)%Veloc(i,0) 
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,0,1) = s_Tdomain%sEdge(ne)%Veloc(i,1) 
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,0,2) = s_Tdomain%sEdge(ne)%Veloc(i,2) 
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(1)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(1) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,0,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,0,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,0,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,0,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,0,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,0,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(2)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(2) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,0,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,0,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,0,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2) 
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,0,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,0,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,0,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(3)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(3) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,0,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,0,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,0,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,0,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,0,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,0,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(4)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(4) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,0,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(5)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(5) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,0,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,0,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(6)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(6) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(0,0,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(0,0,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(0,0,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(0,0,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(7)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(7) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(8)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(8) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,1:s_nglly-2,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(i,0) 
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(i,1) 
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1,s_nglly-1-i,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(i,2) 
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(9)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(9) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(1:s_ngllx-2,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllx-2
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(s_ngllx-1-i,s_nglly-1,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(10)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(10) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,1:s_ngllz-2,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,1:s_ngllz-2,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,1:s_ngllz-2,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_ngllz-2
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1-i,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1-i,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1,s_ngllz-1-i,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif

  ne = s_Tdomain%specel(s_n)%near_edges(11)
  ngll = s_Tdomain%sEdge(ne)%ngll
  if ( s_Tdomain%specel(s_n)%Orient_Edges(11) == 0 ) then
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,0)
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,1)
     s_Tdomain%specel(s_n)%Veloc(0,1:s_nglly-2,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(1:ngll-2,2)
  else
     do i=1,s_nglly-2
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1,0) = s_Tdomain%sEdge(ne)%Veloc(i,0)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1,1) = s_Tdomain%sEdge(ne)%Veloc(i,1)
        s_Tdomain%specel(s_n)%Veloc(0,s_nglly-1-i,s_ngllz-1,2) = s_Tdomain%sEdge(ne)%Veloc(i,2)
     enddo
  endif
end subroutine Veloc_Collection_From_Edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
