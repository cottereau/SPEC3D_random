subroutine define_mechanical_fields_aniso(SubD, &
     ngllx,nglly,ngllz,GlobId,GlobCoord, &
     s_C,yesPML,&
     s_desc,NbElem,NbGlob,nbprocs,n,rg,from_scratch_format,SaveMecaField)


  use ssubdomains

  implicit none

  include 'mpif.h'

  type(subdomain), intent(IN) :: SubD
  integer, intent(IN) :: ngllx, nglly, ngllz, s_desc, nbprocs, n, rg
  integer, dimension(0:nbprocs-1):: NbElem, NbGlob
  integer,dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(IN) ::GlobId
  integer::SZG,SZG23,i0,j0,i1,j1,k1,m, s_code,NbOctInt, NbOctReal
  real,dimension(0:2,0:NbGlob(rg)-1), intent(IN) ::GlobCoord
  logical::yesPML

  integer, dimension(MPI_STATUS_SIZE) :: s_statut
  integer (kind=MPI_OFFSET_KIND) :: PosFile, PosFileI, PosFileR

  integer::desc_chol,gam_ierr 


  integer,parameter::which=1 !indicator of subroutine 'cdfnor'
  real,parameter::mean=0,SD=1 ! ... of Normal germ and 
  real,parameter::X0=0.  ! ... of gaminv
  double precision :: P,Q
  integer::bound,s_status  ! inout of 'gamvinv', 'dpotrf'

  real::delta_diag,Diag,germ
  real,dimension(6,6,0:ngllx-1,0:nglly-1,0:ngllz-1)::c_mean,u_mean,u_fluct
  double precision,dimension(6,6,0:ngllx-1,0:nglly-1,0:ngllz-1),intent(INOUT)::s_C
  double precision,dimension(:,:,:),allocatable::field_m_iem,kappa_eqv,mu_eqv,Kappa_prepost,Mu_prepost!,indice_A
  double precision, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1)::norm_C,indice_A
  real,dimension(21)::Coefs21
  real,dimension(:,:,:,:),allocatable::Field21

!!$  real,dimension(6,6)::u_tempo,c_test
  real,dimension(6,6)::C_temp1,C_temp2,PrePost,Spheric,Deviatoric
  real,dimension(:,:,:,:),allocatable::C2read_write !!!!!!!

  logical::Outer  
  logical,dimension(0:6)::PML_logicals
  !real,dimension(6,6,0:ngllx-1,0:nglly-1,0:ngllz-1)::C_temp3
  logical::Elem_Reload_mat,from_scratch_format,SaveMecaField 


  !print *, 'define_meca::SaveMacaField',SaveMecaField,rg

  if (.not. allocated(C2read_write)) allocate(C2read_write(0:ngllx-1,0:nglly-1,0:ngllz-1,23))

  if (from_scratch_format) then  
     call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,s_code)
     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,s_code)
     SZG = ngllx*nglly*ngllz
     SZG23=23*SZG
     C2read_write(:,:,:,:)=2e38
  
  else

     PML_logicals=.false.
     if (yesPML) then
        PML_logicals(0)=yesPML
        PML_logicals(1)=SubD%Px
        PML_logicals(2)=SubD%Py
        PML_logicals(3)=SubD%Pz
        PML_logicals(4)=SubD%Left
        PML_logicals(5)=SubD%Forward
        PML_logicals(6)=SubD%Down
     endif


     SZG = ngllx*nglly*ngllz
     SZG23=23*SZG



     if (SubD%reload_mat) then
        Elem_Reload_mat = .true.
     else
        Elem_Reload_mat = .false.
     endif


2000 if (SubD%reload_mat .and. Elem_Reload_mat) then

        call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,s_code)
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,s_code)
        PosFile = 2 * nbprocs * NbOctInt ! def_array : mpi_write : Tdomain%n_elem, Tdomain%n_glob_points  
        do i1 = 0, nbprocs-1
           PosFile = PosFile + 3*NbGlob(i1)*NbOctReal ! def_array : mpi_write : Tdomain%GlobCoords
        end do

        PosFileR = PosFile
        PosFileI = PosFile
        do i1 = 0, rg-1 !Previous rg
           PosFileI = PosFileI+NbElem(i1)*NbOctInt*SZG
           PosFileR = PosFileR+NbElem(i1)*NbOctReal*SZG23
        end do

        do i1 = 0, nbprocs-1
           PosFileR = PosFileR + NbElem(i1)*NbOctInt*SZG ! C2read_write "waits for" GlobId
        end do
        PosFileR = PosFileR + n*NbOctReal*SZG23 !C2read_write of previous elts in the same sub-domain
        PosFileI = PosFileI + n*NbOctInt*SZG !GlobId of previous elts in the same sub-domain
        call MPI_FILE_READ_AT(s_desc,PosFileI,GlobId,SZG,MPI_INTEGER,s_statut,s_code)
        call MPI_FILE_READ_AT(s_desc,PosFileR,C2read_write(:,:,:,:),SZG23,MPI_DOUBLE_PRECISION,s_statut,s_code)

        do i0=1,6
           do j0=i0,6
              s_C(i0,j0,:,:,:)=C2read_write(:,:,:,(14-i0)*(i0-1)/2+j0-i0+1)
           enddo
        enddo
!!$print *, minval(abs(s_C))
        if (maxval(s_C) .gt. 2e38) then !The material had not been attributed for this element
           Elem_Reload_mat=.false.
!!$      print *, 'Elem_Reload_false',n, maxval(s_C)
           goto 2000
        endif

     else

        Coefs21=SubD%Coefs21
        !Mean field tensor
        do i0=1,6
           do j0=i0,6
              c_mean(i0,j0,:,:,:)=Coefs21((14-i0)*(i0-1)/2+j0-i0+1)
              if (.not. SubD%SpheDev) u_mean(i0,j0,:,:,:)=c_mean(i0,j0,:,:,:) !Upper stock of c_mean
              if (j0>i0) c_mean(j0,i0,:,:,:)=c_mean(i0,j0,:,:,:)
           enddo
        enddo


        if (.not. SubD%SpheDev) then
           ! Cholesky upper triangle of mean tensor
           do k1 = 0, ngllz-1
              do j1 = 0, nglly-1
                 do i1 = 0, ngllx-1
                    call dpotrf('U',6,u_mean(:,:,i1,j1,k1),6,desc_chol) 
                    !u_mean(:,:,i,j,k) becomes a cholesky's upper triagle of c_mean(:,:,i,j,k) after this operation 
                 enddo
              enddo
           enddo
           if (SubD%random) then      
              allocate(Field21(21,0:ngllx-1,0:nglly-1,0:ngllz-1))
              do m=1,21
                 allocate(field_m_iem(0:ngllx-1,0:nglly-1,0:ngllz-1))
                 call randomisation(SubD%CorNx,SubD%CorNy,SubD%CorNz,&
                      SubD%CorLx,SubD%CorLy,SubD%CorLz,&
                      SubD%ChosenSeed,SubD%Correlation_type,&
                      m,ngllx,nglly,ngllz,GlobCoord,&
                      GlobId,NbGlob,rg,nbprocs,field_m_iem,PML_logicals)
                 Field21(m,:,:,:)=field_m_iem
                 deallocate(field_m_iem)
              enddo
              u_fluct=0.
              germ=0.
              do i0=1,6
                 do j0=i0,6
                    if (j0>i0) then 
                       u_fluct(i0,j0,:,:,:)=sqrt(1./7.)*SubD%delta*Field21(((14-i0)*(i0-1))/2+j0-i0+1,:,:,:)
                    else
                       delta_diag=7/(2*SubD%delta**2)-(2*i0-2)/4
                       do k1 = 0, ngllz-1
                          do j1 = 0, nglly-1
                             do i1 = 0, ngllx-1
                                germ=Field21(((14-i0)*(i0-1))/2+1,i1,j1,k1)
                                select case (SubD%random_type)
                                case("G")
                                   call cdfnor(which,P,Q,germ,mean,SD,s_status,bound)
                                   call gaminv(delta_diag,Diag,X0,P,Q,s_code) 

                                   u_fluct(i0,j0,i1,j1,k1)=sqrt(1./7.)*SubD%delta*sqrt(2*Diag)
                                case("L")
                                   u_fluct(i0,j0,i1,j1,k1)=sqrt(1./7.)*delta_diag**(-0.5)&
                                   *sqrt(2*delta_diag*exp(germ*sqrt(log(delta_diag**(-1)+1))+&
                                        log(1/sqrt(delta_diag**(-1)+1))))                    
                                end select
                             enddo
                          enddo
                       enddo
                    endif
                 enddo
              enddo
              deallocate(Field21)

              do k1 = 0, ngllz-1
                 do j1 = 0, nglly-1
                    do i1 = 0, ngllx-1
                       !The 3 folowing lines are in order to calculate first "transpose(u_mean)*transpose(u_fluct)" then "(.)*u_fluct" and finally "(.)*u_mean"
                       call DGEMM ( 'T', 'T', 6, 6, 6, 1., u_mean(:,:,i1,j1,k1), 6,u_fluct(:,:,i1,j1,k1),6, 0., C_temp1, 6)
                       call DGEMM ( 'N', 'N', 6, 6, 6, 1., C_temp1, 6, u_fluct(:,:,i1,j1,k1),6, 0., C_temp2,6)
                       call DGEMM ( 'N', 'N', 6, 6, 6, 1., C_temp2, 6, u_mean(:,:,i1,j1,k1),6, 0., s_C(:,:,i1,j1,k1),6)
                       !call DGEMM ( 'T', 'N', 6, 6, 6, 1., u_fluct(:,:,i1,j1,k1), 6,u_fluct(:,:,i1,j1,k1),6, 0., s_C(:,:,i1,j1,k1), 6)
                    enddo
                 enddo
              enddo
           else   
              s_C=c_mean
           endif

        else
           !print *,'check1'
           allocate(Kappa_prepost(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate(Mu_prepost(0:ngllx-1,0:nglly-1,0:ngllz-1))
           Mu_prepost=c_mean(4,4,0,0,0)
           Kappa_prepost=2.*c_mean(4,4,0,0,0)/3. + c_mean(1,2,0,0,0)
           Spheric=0.
           Spheric(1:3,1:3)=1.
           Spheric(:,:)=Spheric(:,:)/3.
           Deviatoric(1:6,1:6)=0.
           do i0=1,6
              Deviatoric(i0,i0)=1.
           enddo
           Deviatoric=Deviatoric-Spheric
           !Deviatoric(4:6,4:6)=Deviatoric(4:6,4:6)/sqrt(2.)
           !print *,'check2'
           if (SubD%random) then
              !print *,'check3'
              allocate(Field21(23,0:ngllx-1,0:nglly-1,0:ngllz-1))!N.B. : Field21 stocks 21+2 germs
             !print *,'check4'
              if (.not. (SubD%delta_kernel.lt.1e-4)) then
                 do m=1,21
                    allocate(field_m_iem(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    call randomisation(SubD%CorNx,SubD%CorNy,SubD%CorNz,&
                         SubD%CorLx,SubD%CorLy,SubD%CorLz,&
                         SubD%ChosenSeed,SubD%Correlation_type,&
                         m,ngllx,nglly,ngllz,GlobCoord,&
                         GlobId,NbGlob,rg,nbprocs,field_m_iem,PML_logicals)
                    Field21(m,:,:,:)=field_m_iem
                    deallocate(field_m_iem)
                 enddo
              endif
              if (.not. (SubD%delta.lt.1e-4)) then
                 do m=22,23
                    allocate(field_m_iem(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    call randomisation(SubD%CorNx,SubD%CorNy,SubD%CorNz,&
                         SubD%CorLx,SubD%CorLy,SubD%CorLz,&
                         SubD%ChosenSeed,SubD%Correlation_type,&
                         m,ngllx,nglly,ngllz,GlobCoord,&
                         GlobId,NbGlob,rg,nbprocs,field_m_iem,PML_logicals)
                    Field21(m,:,:,:)=field_m_iem
                    deallocate(field_m_iem)
                 enddo
              endif
              !print *,'check5'
              !Kernel Fluctuation 
              if (.not. (SubD%delta_kernel.lt.1e-4)) then
                 u_fluct=0.
                 germ=0.
                 do i0=1,6
                    do j0=i0,6
                       if (j0>i0) then 
                          u_fluct(i0,j0,:,:,:)=sqrt(1./7.)*SubD%delta_kernel*Field21(((14-i0)*(i0-1))/2+j0-i0+1,:,:,:)
                       else
                          delta_diag=7/(2*SubD%delta_kernel**2)-(2*i0-2)/4
                          do k1 = 0, ngllz-1
                             do j1 = 0, nglly-1
                                do i1 = 0, ngllx-1
                                   germ=Field21(((14-i0)*(i0-1))/2+1,i1,j1,k1)
                                   select case (SubD%random_type_kernel)
                                   case("G")
                                      call cdfnor(which,P,Q,germ,mean,SD,s_status,bound)
                                      call gaminv(delta_diag,Diag,X0,P,Q,s_code) 
                                      u_fluct(i0,j0,i1,j1,k1)=sqrt(1./7.)*SubD%delta_kernel*sqrt(2*Diag)
                                   case("L")
                                      u_fluct(i0,j0,i1,j1,k1)=sqrt(1./7.)*delta_diag**(-0.5)&
                                      *sqrt(2*delta_diag*exp(germ*sqrt(log(delta_diag**(-1)+1))+&
                                           log(1/sqrt(delta_diag**(-1)+1))))                    
                                   end select
                                enddo
                             enddo
                          enddo
                       endif
                    enddo
                 enddo
              endif
              !print *,'check6'
              !Pre-post Fluctuation 
              do k1 = 0, ngllz-1
                 do j1 = 0, nglly-1
                    do i1 = 0, ngllx-1
                       if (.not. (SubD%delta.lt.1e-4)) then


                          delta_diag=1/(SubD%delta**2)
                          !print *,'check61'
                          select case (SubD%random_type)

                          case("G")
                             germ=Field21(22,i1,j1,k1)
                             call cdfnor(which,P,Q,germ,mean,SD,s_status,bound)
                             call gaminv(delta_diag,Diag,X0,P,Q,s_code)
                             Kappa_prepost(i1,j1,k1)=Kappa_prepost(i1,j1,k1)*Diag/delta_diag

                             !print *,'check62'
                             germ=Field21(23,i1,j1,k1)
                             call cdfnor(which,P,Q,germ,mean,SD,s_status,bound)
                             call gaminv(delta_diag,Diag,X0,P,Q,s_code)
                             Mu_prepost(i1,j1,k1)=Mu_prepost(i1,j1,k1)*Diag/delta_diag

                             !print *,'check63'
                          case("L")
                             germ=Field21(22,i1,j1,k1)
                             Kappa_prepost(i1,j1,k1)=Kappa_prepost(i1,j1,k1)*delta_diag&
                             *exp(germ*sqrt(log(delta_diag**(-1)+1))+log(1/sqrt(delta_diag**(-1)+1)))

                             germ=Field21(23,i1,j1,k1)
                             Mu_prepost(i1,j1,k1)=Mu_prepost(i1,j1,k1)*delta_diag&
                             *exp(germ*sqrt(log(delta_diag**(-1)+1))+log(1/sqrt(delta_diag**(-1)+1)))
                          end select
                       endif
                       !print *,'check64'
                       !print *,'check641',sqrt(3*Kappa_prepost(i1,j1,k1)),sqrt(2*Mu_prepost(i1,j1,k1))
                       !print *,'check642',Spheric(1,1),'----------------------'
                       PrePost(:,:)=sqrt(3*Kappa_prepost(i1,j1,k1))*Spheric(:,:)+sqrt(2*Mu_prepost(i1,j1,k1))*Deviatoric(:,:)
                       if (.not. (SubD%delta_kernel.lt.1e-4)) then
                          !print *,'check65'
                          call DGEMM ( 'N', 'T', 6, 6, 6, 1., PrePost(:,:), 6,u_fluct(:,:,i1,j1,k1),6, 0., C_temp1, 6)
                          !print *,'check66'
                          call DGEMM ( 'N', 'N', 6, 6, 6, 1., C_temp1, 6, u_fluct(:,:,i1,j1,k1),6, 0., C_temp2,6)
                          !print *,'check67'
                          call DGEMM ( 'N', 'N', 6, 6, 6, 1., C_temp2, 6, PrePost(:,:),6, 0., s_C(:,:,i1,j1,k1),6)
                          !print *,'check68'
                          !call DGEMM ( 'T', 'N', 6, 6, 6, 1., u_fluct(:,:,i1,j1,k1), 6, u_fluct(:,:,i1,j1,k1),6, 0.,C_temp3,6)
                       else
                          s_C(:,:,i1,j1,k1)=3*Kappa_prepost(i1,j1,k1)*Spheric(:,:)+2*Mu_prepost(i1,j1,k1)*Deviatoric(:,:)
                       endif
                          s_C(4:6,4:6,i1,j1,k1)=s_C(4:6,4:6,i1,j1,k1)/2.
                          s_C(1:3,4:6,i1,j1,k1)=s_C(1:3,4:6,i1,j1,k1)/sqrt(2.)
                          s_C(4:6,1:3,i1,j1,k1)=s_C(4:6,1:3,i1,j1,k1)/sqrt(2.)
                    enddo
                 enddo
              enddo
              deallocate(Field21)
              !print *,'check7'
           else
              s_C=c_mean
           endif
        endif

     endif !test reload




!!!!!!!!!!!!!!!!!!!        KAPPA-MU;Norm; Aniso Indice     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate(kappa_eqv(0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate(mu_eqv(0:ngllx-1,0:nglly-1,0:ngllz-1))

     norm_C=0.
     indice_A=0.
     do k1 = 0, ngllz-1
        do j1 = 0, nglly-1
           do i1 = 0, ngllx-1
              !Condition "outer"

              Outer=PML_logicals(0) .and. &
                   (&
                   ((PML_logicals(1).eqv..true. .and. PML_logicals(4).eqv..true.) .and. (i1<SubD%PML_GLL_out)) .or. &
                   ((PML_logicals(2).eqv..true. .and. PML_logicals(5).eqv..true.) .and. (j1<SubD%PML_GLL_out)) .or. &
                   ((PML_logicals(3).eqv..true. .and. PML_logicals(6).eqv..true.) .and. (k1<SubD%PML_GLL_out)) .or. &
                   ((PML_logicals(1).eqv..true. .and. PML_logicals(4).eqv..false.) .and. (i1>(ngllx-1-SubD%PML_GLL_out))) .or. &
                   ((PML_logicals(2).eqv..true. .and. PML_logicals(5).eqv..false.) .and. (j1>(nglly-1-SubD%PML_GLL_out))) .or. &
                   ((PML_logicals(3).eqv..true. .and. PML_logicals(6).eqv..false.) .and. (k1>(ngllz-1-SubD%PML_GLL_out))) &
                   !!$                ((PML_logicals(1)==.true. .and. PML_logicals(4)==.true.) .and. (float(i1)<(float((ngllx-1))/2))) .or. &
!!$                ((PML_logicals(2)==.true. .and. PML_logicals(5)==.true.) .and. (float(j1)<(float((nglly-1))/2))) .or. &
!!$                ((PML_logicals(3)==.true. .and. PML_logicals(6)==.true.) .and. (float(k1)<(float((ngllz-1))/2))) .or. & 
!!$                ((PML_logicals(1)==.true. .and. PML_logicals(4)==.false.) .and. (float(i1)>(float((ngllx-1))/2))) .or. &
!!$                ((PML_logicals(2)==.true. .and. PML_logicals(5)==.false.) .and. (float(j1)>(float((nglly-1))/2))) .or. &
!!$                ((PML_logicals(3)==.true. .and. PML_logicals(6)==.false.) .and. (float(k1)>(float((ngllz-1))/2))) &
              &)
              kappa_eqv(i1,j1,k1)=(s_C(1,1,i1,j1,k1)+s_C(2,2,i1,j1,k1)+s_C(3,3,i1,j1,k1)&
              +2.*(s_C(1,2,i1,j1,k1)+s_C(1,3,i1,j1,k1)+s_C(2,3,i1,j1,k1)))/9.
              mu_eqv=(s_C(1,1,i1,j1,k1)+s_C(2,2,i1,j1,k1)+s_C(3,3,i1,j1,k1)+3*(s_C(4,4,i1,j1,k1)&
              +s_C(5,5,i1,j1,k1)+s_C(6,6,i1,j1,k1))-(s_C(1,2,i1,j1,k1)+s_C(1,3,i1,j1,k1)+s_C(2,3,i1,j1,k1)))/15.
              !Compute norm(C) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
              do i0=1,6
                 do j0=1,6
                    if (i0>3 .and. j0>3) then
                       norm_C(i1,j1,k1)=norm_C(i1,j1,k1)+4*abs(s_C(i0,j0,i1,j1,k1)-C_mean(i0,j0,i1,j1,k1))**2 !4
                    elseif (i0<4 .and. j0<4) then
                       norm_C(i1,j1,k1)=norm_C(i1,j1,k1)+abs(s_C(i0,j0,i1,j1,k1)-C_mean(i0,j0,i1,j1,k1))**2
                    else
                       norm_C(i1,j1,k1)=norm_C(i1,j1,k1)+2*abs(s_C(i0,j0,i1,j1,k1)-C_mean(i0,j0,i1,j1,k1))**2 !2
                    endif
                 enddo
              enddo
            
              !norm_C(i1,j1,k1)=sqrt(norm_C(i1,j1,k1))

              !C isotropic equivalent & Anisotropy indice
              do i0=1,6
                 do j0=i0,6
                    if (i0==j0) then
                       if (i0>3) then
                          indice_A(i1,j1,k1)=indice_A(i1,j1,k1)+4*(abs(s_C(i0,j0,i1,j1,k1)-mu_eqv(i1,j1,k1)))**2 !4
                          if (Outer) s_C(i0,j0,i1,j1,k1)=mu_eqv(i1,j1,k1)
                       else
                          indice_A(i1,j1,k1)=indice_A(i1,j1,k1)+(abs(s_C(i0,j0,i1,j1,k1)&
                          -4*mu_eqv(i1,j1,k1)/3-kappa_eqv(i1,j1,k1)))**2
                          if (Outer)s_C(i0,j0,i1,j1,k1)=4*mu_eqv(i1,j1,k1)/3+kappa_eqv(i1,j1,k1)
                       endif!i0>3
                    else
                       if (i0>3) then
                          indice_A(i1,j1,k1)=indice_A(i1,j1,k1)+8*(abs(s_C(i0,j0,i1,j1,k1)))**2 !8
                          if (Outer) s_C(i0,j0,i1,j1,k1)=0.
                       else
                          if (j0>3) then
                             indice_A(i1,j1,k1)=indice_A(i1,j1,k1)+4*(abs(s_C(i0,j0,i1,j1,k1)))**2 !4
                             if (Outer) s_C(i0,j0,i1,j1,k1)=0.                     
                          else
                             indice_A(i1,j1,k1)=indice_A(i1,j1,k1)+2*(abs(s_C(i0,j0,i1,j1,k1)&
                             -kappa_eqv(i1,j1,k1)+2*mu_eqv(i1,j1,k1)/3.))**2 !2
                             if (Outer) s_C(i0,j0,i1,j1,k1)=kappa_eqv(i1,j1,k1)-2*mu_eqv(i1,j1,k1)/3.
                          endif
                       endif!i0>3
                       s_C(j0,i0,i1,j1,k1)=s_C(i0,j0,i1,j1,k1)
                    endif !i0=j0
                 enddo
              enddo
              indice_A(i1,j1,k1)=sqrt(indice_A(i1,j1,k1))/norm_C(i1,j1,k1)
              !if (n==1 .and. i1==1 .and. j1==1 .and. k1==1) print *, indice_A(i1,j1,k1)
           enddo
        enddo
     enddo
     deallocate(kappa_eqv,mu_eqv)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!temp 0
     ! allocate(C2read_write(0:ngllx-1,0:nglly-1,0:ngllz-1,23))

     if (SaveMecaField) then   
        do i0=1,6
           do j0=i0,6
              C2read_write(:,:,:,(14-i0)*(i0-1)/2+j0-i0+1)=s_C(i0,j0,:,:,:)
           enddo
        enddo
        !C2read_write(:,:,:,6)=float(rg)
        C2read_write(:,:,:,22)=norm_C(:,:,:)
        C2read_write(:,:,:,23)=indice_A(:,:,:)
     endif
  endif
  !C2read_write(:,:,:,23)=C_temp3(1,1,:,:,:)

  if (SaveMecaField) then 

     call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,s_code)
     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,s_code)
     PosFile = 2 * nbprocs * NbOctInt ! def_array : mpi_write : Tdomain%n_elem, Tdomain%n_glob_points  
     do i1 = 0, nbprocs-1
        PosFile = PosFile + 3*NbGlob(i1)*NbOctReal ! def_array : mpi_write : Tdomain%GlobCoords
     end do

     PosFileR = PosFile
     PosFileI = PosFile
     do i1 = 0, rg-1 !Previous rg
        !print *,NbElem(i1)
        PosFileI = PosFileI+NbElem(i1)*NbOctInt*SZG
        PosFileR = PosFileR+NbElem(i1)*NbOctReal*SZG23
     end do

     do i1 = 0, nbprocs-1
        PosFileR = PosFileR + NbElem(i1)*NbOctInt*SZG ! C2read_write "waits for" GlobId
     end do
     PosFileR = PosFileR + n*NbOctReal*SZG23 !C2read_write of previous elts in the same sub-domain
     PosFileI = PosFileI + n*NbOctInt*SZG !GlobId of previous elts in the same sub-domain
     call MPI_FILE_WRITE_AT(s_desc,PosFileI,GlobId,SZG,MPI_INTEGER,s_statut,s_code)
     call MPI_FILE_WRITE_AT(s_desc,PosFileR,C2read_write(:,:,:,:),SZG23,MPI_DOUBLE_PRECISION,s_statut,s_code)
     !print *, 'PosFileI,PosFileR', PosFileI, PosFileR
  endif
  deallocate(C2read_write)

end subroutine define_mechanical_fields_aniso




!----------------------------------------------------------

!SUBROUTINE generating the stochastic field of elasticity tensor's germ 
subroutine randomisation(CorNx,CorNy,CorNz,CorLx,CorLy,CorLz,ChosenSeed,&
     Correlation_type,s_m,s_ngllx,s_nglly,s_ngllz,s_GlobCoord,&
     s_GlobId,s_NbGlob,s_rg,s_nbprocs,s_RandParam,s_PML_logicals)
  use pig  
  implicit none
  integer::s_m,s_ngllx,s_nglly,s_ngllz,ic, jc, kc, s_code, i, j, k, CorNx, CorNy, CorNz,ChosenSeed,s_rg,s_nbprocs
!!$  integer, dimension(0:2*CorNx*CorNy*CorNz-1) :: SEED
  integer,dimension(:),allocatable::SEED
  integer,dimension(0:s_ngllx-1,0:s_nglly-1,0:s_ngllz-1), intent(IN) ::s_GlobId
  integer, dimension(0:s_nbprocs-1)::s_NbGlob


  real,dimension(0:2,0:s_NbGlob(s_rg)-1), intent(IN) ::s_GlobCoord
  real,dimension(0:s_ngllx-1,0:s_nglly-1,0:s_ngllz-1),intent(INOUT)::s_RandParam
!!$  real,dimension(0:s_ngllx-1,0:s_nglly-1,0:s_ngllz-1,0:1)::RandParamUnif
  real,dimension(:,:,:,:),allocatable::RandParamUnif
  real,dimension(0:CorNx*CorNy*CorNz-1)::tripl0,tripl1,tripl2,dot_prod,S
  real,dimension(0:CorNx*CorNy*CorNz-1,0:2)::K_support !support of S  
  real :: dKx, dKy, dKz, Xc, Yc, Zc, CorLx, CorLy, CorLz
  character(len=1)::Correlation_type

  logical,dimension(0:6)::s_PML_logicals

  if (.not. allocated (RandParamUnif)) allocate(RandParamUnif(0:CorNx-1,0:CorNy-1,0:CorNz-1,0:1))
  if (.not. allocated (SEED)) allocate(SEED(0:2*CorNx*CorNy*CorNz-1))

  dKx = 2*pi/CorLx/CorNx
  dKy = 2*pi/CorLy/CorNy
  dKz = 2*pi/CorLz/CorNz

  do i = 0, CorNx-1
     do j = 0, CorNy-1
        do k = 0, CorNz-1
           K_support(i*CorNy*CorNz+j*CorNz+k,0)=-pi/CorLx+(i+1)*dKx 
           K_support(i*CorNy*CorNz+j*CorNz+k,1)=-pi/CorLy+(j+1)*dKy
           K_support(i*CorNy*CorNz+j*CorNz+k,2)=-pi/CorLz+(k+1)*dKz
           do kc = 1, 2
              SEED(i+CorNx*(j+CorNy*(k+CorNz*(kc-1))))= i+CorNx*(j+CorNy*(k+CorNz*(kc-1)))+ &
                   ChosenSeed*23+s_m 
           enddo
        enddo
     enddo
  enddo
  S=K_support(:,0)
  call tripuls(S,size(S),pi/CorLx,tripl0,Correlation_type)


  S=K_support(:,1)
  call tripuls(S,size(S),pi/CorLy,tripl1,Correlation_type)


  S=K_support(:,2)
  call tripuls(S,size(S),pi/CorLz,tripl2,Correlation_type)


  call arrmul(tripl0,tripl1,tripl2,size(S),dot_prod)



  S=(CorLx/pi)*(CorLy/pi)*(CorLz/pi)*dot_prod


  call random_seed(PUT=SEED(0:2*CorNx*CorNy*CorNz-1))
  call random_number(RandParamUnif)

  RandParamUnif(:,:,:,0)=sqrt(-log(RandParamUnif(:,:,:,0))) !amplitude
  RandParamUnif(:,:,:,1)=RandParamUnif(:,:,:,1)*2*pi !phase

  do k = 0, s_ngllz-1
     do j = 0, s_nglly-1
        do i = 0, s_ngllx-1

           s_RandParam(i,j,k)= 0
           if (.not. s_PML_logicals(0)) then
              Xc = s_GlobCoord(0,s_GlobId(i,j,k))
              Yc = s_GlobCoord(1,s_GlobId(i,j,k))
              Zc = s_GlobCoord(2,s_GlobId(i,j,k))
           else
              if (((s_PML_logicals(1)).and.(.not.s_PML_logicals(2))).and.(.not.s_PML_logicals(3))) then
                 if (s_PML_logicals(4)) then
                    Xc = s_GlobCoord(0,s_GlobId(s_ngllx-1,j,k))
                    Yc = s_GlobCoord(1,s_GlobId(s_ngllx-1,j,k))
                    Zc = s_GlobCoord(2,s_GlobId(s_ngllx-1,j,k))
                 else
                    Xc = s_GlobCoord(0,s_GlobId(0,j,k))
                    Yc = s_GlobCoord(1,s_GlobId(0,j,k))
                    Zc = s_GlobCoord(2,s_GlobId(0,j,k))
                 endif
              endif

              if (((s_PML_logicals(2)).and.(.not.s_PML_logicals(1))).and.(.not.s_PML_logicals(3))) then
                 if (s_PML_logicals(5)) then
                    Xc = s_GlobCoord(0,s_GlobId(i,s_nglly-1,k))
                    Yc = s_GlobCoord(1,s_GlobId(i,s_nglly-1,k))
                    Zc = s_GlobCoord(2,s_GlobId(i,s_nglly-1,k))
                 else
                    Xc = s_GlobCoord(0,s_GlobId(i,0,k))
                    Yc = s_GlobCoord(1,s_GlobId(i,0,k))
                    Zc = s_GlobCoord(2,s_GlobId(i,0,k))
                 endif
              endif

              if (((s_PML_logicals(3)).and.(.not.s_PML_logicals(1))).and.(.not.s_PML_logicals(2))) then
                 if (s_PML_logicals(6)) then
                    Xc = s_GlobCoord(0,s_GlobId(i,j,s_ngllz-1))
                    Yc = s_GlobCoord(1,s_GlobId(i,j,s_ngllz-1))
                    Zc = s_GlobCoord(2,s_GlobId(i,j,s_ngllz-1))
                 else
                    Xc = s_GlobCoord(0,s_GlobId(i,j,0))
                    Yc = s_GlobCoord(1,s_GlobId(i,j,0))
                    Zc = s_GlobCoord(2,s_GlobId(i,j,0))
                 endif
              endif

           endif

           do kc = 0, CorNz-1
              do jc = 0, CorNy-1
                 do ic = 0, CorNx-1
                    !germ construction
                    s_RandParam(i,j,k) = s_RandParam(i,j,k) + &
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

end subroutine randomisation
!-------------------------------------------------------------------------------------------

!SUBROUTINE get symmetric-unity amplitude tripuls value
subroutine tripuls(x,size_x,half_width,y,Corr_type)
  implicit none
  real::half_width
  integer::i,size_x
  real,dimension(size_x)::x
  real,dimension(size_x),intent(OUT)::y
  character(len=1) :: Corr_type !Correlation_type
  !real,dimension(:),allocatable::tripuls
  !real,dimension(size(x))::tripuls

  select case (Corr_type)
  case ("S") ! sinc^2
     do i=1,size(x)
        if ((abs(x(i)) .GT. half_width) .or. (abs(x(i)) == half_width) ) then
           y(i)=0
        else 
           y(i)=1-(1/half_width)*abs(x(i))
        endif
     enddo
  case ("E") ! Exponential
     do i=1,size(x)
        y(i)=exp(-2*abs(x(i))/half_width/sqrt(3.14159))
     enddo
  end select
end subroutine tripuls
!----------------------------------------------------------

!SUBROUTINE compute vector multiplication of type A.*B.*C in matlab
subroutine arrmul(x1,x2,x3,size_x,z)
  implicit none
  integer::j,size_x
  real,dimension(size_x),intent(IN)::x1,x2,x3
  real,dimension(size_x),intent(OUT)::z
!!$!real,dimension(:),allocatable::arrmul
!!$!real,dimension(size(x))::arrmul

  !if (.not. allocated (y)) allocate(y(0:size(x1)-1))
  do j=1,size(x1)
     z(j)=x1(j)*x2(j)*x3(j)
  enddo
end subroutine arrmul
!----------------------------------------------------------------






!$Id: define_mechanical_fields.f90,v 1.3 2007/08/30 15:40:55 taquanganh Exp $!  


