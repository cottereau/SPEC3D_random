subroutine define_mechanical_fields_iso(SubD, &
     ngllx,nglly,ngllz,GlobId,GlobCoord, &
     RLam,RMu,lPML, &
     s_desc,NbElem,NbGlob,nbprocs,n,rg)

  use pig
  use ssubdomains

  implicit none

  include 'mpif.h'

  type(subdomain), intent(IN) :: SubD
  integer, intent(IN) :: ngllx, nglly, ngllz, s_desc, nbprocs, n, rg
  integer, dimension(0:nbprocs-1), intent(IN) :: NbElem, NbGlob
  integer,dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(IN) ::GlobId
  real,dimension(0:2,0:NbGlob(rg)-1), intent(IN) ::GlobCoord
  real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1),intent(OUT) :: RLam,RMu
  logical,intent(IN):: lPML

  real,dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: RandParam,RandParam_trans
!!$  real, dimension(0:SubD%CorNx-1,0:SubD%CorNy-1,0:SubD%CorNz-1,1:4) :: RandParamUnif
  real, dimension(:,:,:,:),allocatable :: RandParamUnif
  real :: dKx, dKy, dKz, Ampl, Xc, Yc, Zc, CorLx, CorLy, CorLz, delta, MeanParam
  integer :: SZG, ic, jc, kc, s_code, i, j, k, CorNx, CorNy, CorNz, NbOctInt, NbOctReal
!!$  integer, dimension(0:4*SubD%CorNx*SubD%CorNy*SubD%CorNz-1) :: SEED
  integer, dimension(:),allocatable :: SEED
  integer, dimension(MPI_STATUS_SIZE) :: s_statut
  integer (kind=MPI_OFFSET_KIND) :: PosFile, PosFileI, PosFileR

!!$  real,dimension(:),allocatable::tripl0,tripl1,tripl2,dot_prod,S !S:Elementery Spectral density Vector (reshape)
!!$  real,dimension(:,:),allocatable::K_support !support of S  
  real,dimension(0:SubD%CorNx*SubD%CorNy*SubD%CorNz-1)::tripl0,tripl1,tripl2,dot_prod,S
  real,dimension(0:SubD%CorNx*SubD%CorNy*SubD%CorNz-1,0:2)::K_support !support of S  
  integer,parameter::which=1 !indicator of subroutine 'cdfnor'
  real,parameter::mean=0,SD=1 ! ... of Normal germ and 
  real,parameter::X0=0  ! ... of gaminv
  double precision :: P,Q
  integer::bound,s_status  ! inout of 'gamvinv'

  Rlam = SubD%DLambda
  Rmu  = SubD%DMu
  if (SubD%random) then      

     CorNx = SubD%CorNx
     CorNy = SubD%CorNy
     CorNz = SubD%CorNz
     CorLx = SubD%CorLx
     CorLy = SubD%CorLy
     CorLz = SubD%CorLz
     delta = SubD%delta

     select case (SubD%random_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case ("N") ! Normal law for local value
        if (.not. allocated (RandParamUnif)) allocate(RandParamUnif(0:CorNx-1,0:CorNy-1,0:CorNz-1,1:4))
        if (.not. allocated (SEED)) allocate(SEED(0:4*CorNx*CorNy*CorNz-1))
        dKx = 5/CorLx/CorNx
        dKy = 5/CorLy/CorNy
        dKz = 5/CorLz/CorNz
        Ampl=delta*sqrt(CorLx*CorLy*CorLz*dKx*dKy*dKz/pi)

        do i = 0, CorNx-1
           do j = 0, CorNy-1
              do k = 0, CorNz-1
                 do kc = 1, 4
                    SEED(i+CorNx*(j+CorNy*(k+CorNz*(kc-1))))= i+CorNx*(j+CorNy*(k+CorNz*(kc-1)))+ &
                         SubD%ChosenSeed
                 end do
              end do
           end do
        end do
        call random_seed(PUT=SEED(0:4*CorNx*CorNy*CorNz-1))

        call random_number(RandParamUnif)
        RandParamUnif=RandParamUnif*2*pi

        do k = 0, ngllz-1
           do j = 0, nglly-1
              do i = 0, ngllx-1
                 RandParam(i,j,k)= 0
                 Xc = GlobCoord(0,GlobId(i,j,k))
                 Yc = GlobCoord(1,GlobId(i,j,k))
                 Zc = GlobCoord(2,GlobId(i,j,k))
                 do kc = 0, CorNz-1
                    do jc = 0, CorNy-1
                       do ic = 0, CorNx-1
                          RandParam(i,j,k) = RandParam(i,j,k) + &
                               sqrt(exp(-(ic*dKx*CorLx/2)**2-(jc*dKy*CorLy/2)**2-(kc*dKz*CorLz/2)**2))*( &
                               cos(ic*dKx*Xc+jc*dKy*Yc+kc*dKz*Zc+RandParamUnif(ic,jc,kc,1)) + &
                               cos(ic*dKx*Xc+jc*dKy*Yc-kc*dKz*Zc+RandParamUnif(ic,jc,kc,2)) + &
                               cos(ic*dKx*Xc-jc*dKy*Yc+kc*dKz*Zc+RandParamUnif(ic,jc,kc,3)) + &
                               cos(ic*dKx*Xc-jc*dKy*Yc-kc*dKz*Zc+RandParamUnif(ic,jc,kc,4)))
!!$                           Code par Regis 05/2007                        
!!$                           M. Shinozuka and G. Deodatis. Simulation of multidimensional gaussian 
!!$                           stochastic fields by spectral representation. 49(1):29-53, 1996.
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
        deallocate(RandParamUnif,SEED)
        RandParam=1+RandParam*Ampl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case ("G") ! Gamma law for local value of Rmu
        !delta = Tdomain%sSubdomain(Tdomain%specel(n)%mat_index)%delta !coeficient of variation
        dKx = 2*pi/CorLx/CorNx
        dKy = 2*pi/CorLy/CorNy
        dKz = 2*pi/CorLz/CorNz
!!$        SZG=2*CorNx*CorNy*CorNz !2 :: number of random var in each element of the trigonometric base :: phase & amplitude

        if (.not. allocated (RandParamUnif)) allocate(RandParamUnif(0:CorNx-1,0:CorNy-1,0:CorNz-1,0:1))
        if (.not. allocated (SEED)) allocate(SEED(0:2*CorNx*CorNy*CorNz-1))

        do i = 0, CorNx-1
           do j = 0, CorNy-1
              do k = 0, CorNz-1
                 K_support(i*CorNy*CorNz+j*CorNz+k,0)=-pi/CorLx+(i+1)*dKx 
                 K_support(i*CorNy*CorNz+j*CorNz+k,1)=-pi/CorLy+(j+1)*dKy
                 K_support(i*CorNy*CorNz+j*CorNz+k,2)=-pi/CorLz+(k+1)*dKz
                 do kc = 1, 2
                    SEED(i+CorNx*(j+CorNy*(k+CorNz*(kc-1))))= i+CorNx*(j+CorNy*(k+CorNz*(kc-1)))+ &
                         SubD%ChosenSeed 
                 enddo
              enddo
           enddo
        enddo
        S=K_support(:,0)
        call tripuls(S,size(S),pi/CorLx,tripl0,SubD%Correlation_type)


        S=K_support(:,1)
        call tripuls(S,size(S),pi/CorLy,tripl1,SubD%Correlation_type)


        S=K_support(:,2)
        call tripuls(S,size(S),pi/CorLz,tripl2,SubD%Correlation_type)


        call arrmul(tripl0,tripl1,tripl2,size(S),dot_prod)



        S=(CorLx/pi)*(CorLy/pi)*(CorLz/pi)*dot_prod



        call random_seed(PUT=SEED(0:2*CorNx*CorNy*CorNz-1))
        call random_number(RandParamUnif)

        RandParamUnif(:,:,:,0)=sqrt(-log(RandParamUnif(:,:,:,0))) !amplitude
        RandParamUnif(:,:,:,1)=RandParamUnif(:,:,:,1)*2*pi !phase


        do k = 0, ngllz-1
           do j = 0, nglly-1
              do i = 0, ngllx-1

                 RandParam(i,j,k)= 0
                 Xc = GlobCoord(0,GlobId(i,j,k))
                 Yc = GlobCoord(1,GlobId(i,j,k))
                 Zc = GlobCoord(2,GlobId(i,j,k))
                 do kc = 0, CorNz-1
                    do jc = 0, CorNy-1
                       do ic = 0, CorNx-1
                          !germ construction
                          if (ic==15 .and. jc==0 .and. kc==0 .and. n==0) then
                             print *,sqrt(2*dKx*dKy*dKz), sqrt(S(ic*CorNy*CorNz+jc*CorNz+kc)),S(ic*CorNy*CorNz+jc*CorNz+kc)  
                             print *, K_support(ic*CorNy*CorNz+jc*CorNz+kc,0),Xc 
                             print *, K_support(ic*CorNy*CorNz+jc*CorNz+kc,1),Yc
                             print *, K_support(ic*CorNy*CorNz+jc*CorNz+kc,2),Zc
                          endif
                          RandParam(i,j,k) = RandParam(i,j,k) + &
                               sqrt(2*dKx*dKy*dKz)*sqrt(S(ic*CorNy*CorNz+jc*CorNz+kc))*&
                               RandParamUnif(ic,jc,kc,0)*&
                               cos(RandParamUnif(ic,jc,kc,1)+K_support(ic*CorNy*CorNz+jc*CorNz+kc,0)*Xc+&
                               K_support(ic*CorNy*CorNz+jc*CorNz+kc,1)*Yc+&
                               K_support(ic*CorNy*CorNz+jc*CorNz+kc,2)*Zc)
                          !print *,RandParam(i,j,k), ic, jc, kc
                          !gamma transformation                         
                          call cdfnor(which,P,Q,RandParam(i,j,k),mean,SD,s_status,bound)
                          call gaminv(delta**(-2),RandParam_trans(i,j,k),X0,P,Q,s_code)
                          !print *,RandParam_trans(i,j,k), s_code         
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        RandParam=RandParam_trans*delta**2
        deallocate(RandParamUnif,SEED)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     case ("L") ! Log-Normal law for local value of Rmu

        dKx = 2*pi/CorLx/CorNx
        dKy = 2*pi/CorLy/CorNy
        dKz = 2*pi/CorLz/CorNz


        if (.not. allocated (RandParamUnif)) allocate(RandParamUnif(0:CorNx-1,0:CorNy-1,0:CorNz-1,0:1))

        if (.not. allocated (SEED)) allocate(SEED(0:2*CorNx*CorNy*CorNz-1))
        do i = 0, CorNx-1
           do j = 0, CorNy-1
              do k = 0, CorNz-1
                 K_support(i*CorNy*CorNz+j*CorNz+k,0)=-pi/CorLx+(i+1)*dKx 
                 K_support(i*CorNy*CorNz+j*CorNz+k,1)=-pi/CorLy+(j+1)*dKy
                 K_support(i*CorNy*CorNz+j*CorNz+k,2)=-pi/CorLz+(k+1)*dKz
                 do kc = 1, 2
                    SEED(i+CorNx*(j+CorNy*(k+CorNz*(kc-1))))= i+CorNx*(j+CorNy*(k+CorNz*(kc-1)))+ &
                         SubD%ChosenSeed 
                 enddo
              enddo
           enddo
        enddo
        S=K_support(:,0)
        call tripuls(S,size(S),pi/CorLx,tripl0,SubD%Correlation_type)


        S=K_support(:,1)
        call tripuls(S,size(S),pi/CorLy,tripl1,SubD%Correlation_type)


        S=K_support(:,2)
        call tripuls(S,size(S),pi/CorLz,tripl2,SubD%Correlation_type)


        call arrmul(tripl0,tripl1,tripl2,size(S),dot_prod)



        S=(CorLx/pi)*(CorLy/pi)*(CorLz/pi)*dot_prod


        call random_seed(PUT=SEED(0:2*CorNx*CorNy*CorNz-1))
        call random_number(RandParamUnif)

        RandParamUnif(:,:,:,0)=sqrt(-log(RandParamUnif(:,:,:,0))) !amplitude
        RandParamUnif(:,:,:,1)=RandParamUnif(:,:,:,1)*2*pi !phase

        do k = 0, ngllz-1
           do j = 0, nglly-1
              do i = 0, ngllx-1

                 RandParam(i,j,k)= 0
                 Xc = GlobCoord(0,GlobId(i,j,k))
                 Yc = GlobCoord(1,GlobId(i,j,k))
                 Zc = GlobCoord(2,GlobId(i,j,k))
                 do kc = 0, CorNz-1
                    do jc = 0, CorNy-1
                       do ic = 0, CorNx-1
                          !germ construction
                          RandParam(i,j,k) = RandParam(i,j,k) + &
                               sqrt(2*dKx*dKy*dKz)*sqrt(S(ic*CorNy*CorNz+jc*CorNz+kc))*&
                               RandParamUnif(ic,jc,kc,0)*&
                               cos(RandParamUnif(ic,jc,kc,1)+K_support(ic*CorNy*CorNz+jc*CorNz+kc,0)*Xc+&
                               K_support(ic*CorNy*CorNz+jc*CorNz+kc,1)*Yc+&
                               K_support(ic*CorNy*CorNz+jc*CorNz+kc,2)*Zc)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        deallocate(RandParamUnif,SEED)
        !Log-Normal transformation 
        RandParam=exp(RandParam*sqrt(log(delta**2+1))+log(1/sqrt(delta**2+1))) 

     end select


     if (.not.lPML) then
        Rlam = RLam * RandParam
        Rmu = RMu * RandParam
     elseif (((SubD%Px).and.(.not.SubD%Py)).and.(.not.SubD%Pz)) then
        MeanParam=0
        if (SubD%Left) then
           do k = 0, ngllz-1
              do j = 0, nglly-1
                 MeanParam = MeanParam + RandParam(ngllx-1,j,k)
              enddo
           enddo
        else
           do k = 0, ngllz-1
              do j = 0, nglly-1
                 MeanParam = MeanParam + RandParam(0,j,k)
              enddo
           enddo
        endif
        Rlam = RLam * MeanParam/(nglly*ngllz)
        Rmu  = RMu * MeanParam/(nglly*ngllz)
     elseif (((SubD%Py).and.(.not.SubD%Px)).and.(.not.SubD%Pz)) then
        MeanParam=0
        if (SubD%Forward) then
           do k = 0, ngllz-1
              do i = 0, ngllx-1
                 MeanParam = MeanParam + RandParam(i,nglly-1,k)
              enddo
           enddo
        else
           do k = 0, ngllz-1
              do i = 0, ngllx-1
                 MeanParam = MeanParam + RandParam(i,0,k)
              enddo
           enddo
        endif
        Rlam = RLam * MeanParam/(ngllx*ngllz)
        Rmu  = RMu * MeanParam/(ngllx*ngllz)
     elseif (((SubD%Pz).and.(.not.SubD%Px)).and.(.not.SubD%Py)) then
        MeanParam=0
        if (SubD%Down) then
           do j = 0, nglly-1
              do i = 0, ngllx-1
                 MeanParam = MeanParam + RandParam(i,j,ngllz-1)
              enddo
           enddo
        else
           do j = 0, nglly-1
              do i = 0, ngllx-1
                 MeanParam = MeanParam + RandParam(i,j,0)
              enddo
           enddo
        endif
        Rlam = RLam * MeanParam/(nglly*ngllx)
        Rmu  = RMu * MeanParam/(nglly*ngllx)
     endif

  endif

  if (.not. (s_desc==-1)) then
     SZG = ngllx*nglly*ngllz
     call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,s_code)
     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,s_code)
     PosFile = 2 * nbprocs * NbOctInt
     do i = 0, nbprocs-1
        PosFile = PosFile + 3*NbGlob(i)*NbOctReal
     end do
     PosFileR = PosFile
     PosFileI = PosFile
     do i = 0, rg-1
        PosFileI = PosFileI+NbElem(i)*NbOctInt*SZG
        PosFileR = PosFileR+NbElem(i)*NbOctReal*SZG
     end do
     do i = 0, nbprocs-1
        PosFileR = PosFileR + NbElem(i)*NbOctInt*SZG
     end do
     PosFileR = PosFileR + n*NbOctReal*SZG
     PosFileI = PosFileI + n*NbOctInt*SZG
     call MPI_FILE_WRITE_AT(s_desc,PosFileI,GlobId,SZG,MPI_INTEGER,s_statut,s_code)
     call MPI_FILE_WRITE_AT(s_desc,PosFileR,Rmu,SZG,MPI_DOUBLE_PRECISION,s_statut,s_code) 
  endif
end subroutine define_mechanical_fields_iso



!!$
!!$
!!$!SUBROUTINE get symmetric-unity amplitude tripuls value
!!$subroutine tripuls(x,size_x,half_width,y,Corr_type)
!!$  implicit none
!!$  real::half_width
!!$  integer::i,size_x
!!$  real,dimension(size_x)::x
!!$  real,dimension(size_x),intent(OUT)::y
!!$  character(len=1) :: Corr_type !Correlation_type
!!$  !real,dimension(:),allocatable::tripuls
!!$  !real,dimension(size(x))::tripuls
!!$
!!$  select case (Corr_type)
!!$  case ("S") ! sinc^2
!!$     do i=0,size(x)-1
!!$        if ((abs(x(i)) .GT. half_width) .or. (abs(x(i)) == half_width) ) then
!!$           y(i)=0
!!$        else 
!!$           y(i)=1-(1/half_width)*abs(x(i))
!!$        endif
!!$     enddo
!!$  case ("E") ! Exponential
!!$     do i=0,size(x)-1
!!$        y(i)=exp(-2*abs(x(i))/half_width/sqrt(3.14159))
!!$     enddo
!!$  end select
!!$end subroutine tripuls
!!$!----------------------------------------------------------
!!$
!!$!SUBROUTINE compute vector multiplication of type A.*B.*C in matlab
!!$subroutine arrmul(x1,x2,x3,size_x,z)
!!$  implicit none
!!$  integer::j,size_x
!!$  real,dimension(size_x),intent(IN)::x1,x2,x3
!!$  real,dimension(size_x),intent(OUT)::z
!real,dimension(:),allocatable::arrmul
!real,dimension(size(x))::arrmul
!!$
!!$  !if (.not. allocated (y)) allocate(y(0:size(x1)-1))
!!$  do j=0,size(x1)-1
!!$     z(j)=x1(j)*x2(j)*x3(j)
!!$  enddo
!!$end subroutine arrmul
!!$!----------------------------------------------------------------
!!$
!!$




!$Id: define_mechanical_fields.f90,v 1.3 2007/08/30 15:40:55 taquanganh Exp $!  


