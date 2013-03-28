module selement

  implicit none

  ! Modified by Gaetano Festa 23/02/2005

  type :: element

     integer :: mat_index, ngllx, nglly, ngllz
     integer, dimension (:), pointer :: Control_nodes
      integer, dimension (0:5) :: Near_Faces, Orient_Faces
     integer, dimension (0:11) :: Near_Edges, Orient_Edges
     integer, dimension (0:7) :: Near_Vertices
     integer, dimension (:,:,:), pointer :: Iglobnum

     real, dimension (:,:,:), pointer :: Jacob, Density, Lambda, Mu, MassMat,MaxAbsVeloc
     real, dimension(:,:,:,:), pointer :: ACoeff, Forces,Veloc,Displ,Accel,V0, Diagonal_Stress, Residual_Stress,stress_sigma,TravelTime
     real, dimension(:,:,:,:,:), pointer :: InvGrad,c
     logical, dimension(:,:,:,:),pointer :: TravelTimeFound
!!$     real, dimension(:,:,:,:), allocatable :: D0
!!$     integer,dimension(:),pointer::signus

     ! Anisotropy allocation
     logical :: aniso 

     ! PML allocation 
     logical :: PML,FPML
     real, dimension (:,:,:,:), pointer :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3, & 
                                                             Diagonal_Stress4a,Diagonal_Stress4b, &
                                                             Diagonal_Stress5a,Diagonal_Stress5b, &
                                                             Diagonal_Stress6a,Diagonal_Stress6b
     real, dimension (:,:,:,:), pointer :: Residual_Stress1, Residual_Stress2,Residual_stress3, &
                                                             Residual_Stress4a,Residual_stress4b, &
                                                             Residual_Stress5a,Residual_stress5b, &
                                                             Residual_Stress6a,Residual_stress6b
     real, dimension (:,:,:,:), pointer :: DumpSx,DumpSy,DumpSz
     real, dimension (:,:,:,:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3
     real, dimension (:,:,:,:), pointer :: DumpVx,DumpVy,DumpVz, DumpMass

     real, dimension  (:,:,:,:), pointer :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
     real, dimension  (:,:,:,:), pointer :: I_Residual_Stress1, I_Residual_Stress2
     real, dimension  (:,:,:,:), pointer :: Iveloc1, Iveloc2, Iveloc3
     real, dimension (:,:,:), pointer :: Isx, Isy, Isz, Ivx, Ivy, Ivz

     !logical, reload_mat
  end type element

contains 

  ! ############################################################
  subroutine Prediction_Elem_Veloc (Elem, alpha, bega, gam1, dt)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: alpha, bega, gam1, dt

    integer :: ngllx, nglly, ngllz


    ngllx = Elem%ngllx ; nglly =Elem%nglly;  ngllz = Elem%ngllz

    Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) &
    + dt * Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) + dt**2 * (0.5 - bega) * Elem%Accel
    !U_{n+1}=U_{n} + dt*V_{n} + dt^2(0.5-beta/gamma)*A_{n}  
    Elem%V0 = Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)
    Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = alpha * Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) + (1-alpha) * Elem%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) 
    !U_{n+alpha}=alpha* U_{n+1}+(1-alpha)*U_{n}
!!$    if (allocated(Elem%D0)) Elem%D0=Elem%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)
    return
  end subroutine Prediction_Elem_Veloc

  ! ###########################################################
  subroutine Correction_Elem_Veloc (Elem, bega, gam1, dt)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: bega, gam1, dt

    integer :: ngllx, nglly, ngllz, i


    ngllx = Elem%ngllx ;  nglly =Elem%nglly ;  ngllz = Elem%ngllz
    do i = 0,2
       Elem%Forces(1:ngllx-2,1:nglly-2, 1:ngllz-2,i) = Elem%MassMat(:,:,:)  * Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
    enddo
    Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)  = Elem%v0(:,:,:,:)+ dt * Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,:)
    Elem%Accel  = Elem%Accel + gam1 /dt * (Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)-Elem%V0)
    Elem%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) &
    + bega *dt * (Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)+Elem%V0)

    return
  end subroutine Correction_Elem_Veloc

  ! ###########################################################
  subroutine  compute_InternalForces_Elem (Elem, hprimex, htprimex, hprimey, htprimey, hprimez, htprimez)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex, htprimex
    real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
    real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez

    integer :: n_z, m1, m2, m3

    real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
    real, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
    real, dimension ( 0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
         dUy_dxi, dUy_deta, dUy_dzeta, &
         dUz_dxi, dUz_deta, dUz_dzeta, &
         t1, s0, Uxloc, Uyloc, Uzloc

    m1 = Elem%ngllx; m2 = Elem%nglly; m3= Elem%ngllz

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,0) ,m1, 0., dUx_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,0), m1, hprimey ,m2, 0., dUx_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0., dUx_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,1) ,m1, 0., dUy_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,1), m1, hprimey ,m2, 0., dUy_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0., dUy_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,2) ,m1, 0., dUz_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,2), m1, hprimey ,m2, 0., dUz_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0., dUz_dzeta, m1*m2 )

    t1 = Elem%Acoeff(:,:,:,0)*dUx_dxi + Elem%Acoeff(:,:,:,1)*dUx_deta + Elem%Acoeff(:,:,:,2)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,3)*dUy_dxi + Elem%Acoeff(:,:,:,4)*dUy_deta + Elem%Acoeff(:,:,:,5)*dUy_dzeta +  &
         Elem%Acoeff(:,:,:,6)*dUz_dxi + Elem%Acoeff(:,:,:,7)*dUz_deta + Elem%Acoeff(:,:,:,8)*dUz_dzeta 

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,t1(0,0,0) ,m1, 0., Uxloc, m1 )

    t1 = Elem%Acoeff(:,:,:,1)*dUx_dxi + Elem%Acoeff(:,:,:,9)*dUx_deta + Elem%Acoeff(:,:,:,10)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,11)*dUy_dxi + Elem%Acoeff(:,:,:,12)*dUy_deta + Elem%Acoeff(:,:,:,13)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,14)*dUz_dxi + Elem%Acoeff(:,:,:,15)*dUz_deta + Elem%Acoeff(:,:,:,16)*dUz_dzeta

    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1.,t1(0,0,n_z), m1, htprimey ,m2, 0., s0(0,0,n_z),m1 )
    enddo
    Uxloc = s0 + Uxloc

    t1 = Elem%Acoeff(:,:,:,2)*dUx_dxi + Elem%Acoeff(:,:,:,10)*dUx_deta + Elem%Acoeff(:,:,:,17)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,18)*dUy_dxi + Elem%Acoeff(:,:,:,19)*dUy_deta + Elem%Acoeff(:,:,:,20)*dUy_dzeta +&
         Elem%Acoeff(:,:,:,21)*dUz_dxi + Elem%ACoeff(:,:,:,22)*dUz_deta + Elem%Acoeff(:,:,:,23)*dUz_dzeta

    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., t1(0,0,0), m1*m2, htprimez ,m3, 0., s0, m1*m2 )
    Uxloc = s0 + Uxloc

    t1 = Elem%Acoeff(:,:,:,3)*dUx_dxi + Elem%Acoeff(:,:,:,11)*dUx_deta + Elem%Acoeff(:,:,:,18)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,24)*dUy_dxi + Elem%Acoeff(:,:,:,25)*dUy_deta + Elem%Acoeff(:,:,:,26)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,27)*dUz_dxi + Elem%Acoeff(:,:,:,28)*dUz_deta + Elem%Acoeff(:,:,:,29)*dUz_dzeta

    call DGEMM ('N', 'N', m1, m2*m3, m1, 1., hprimex, m1, t1(0,0,0), m1, 0., Uyloc, m1)

    t1 = Elem%Acoeff(:,:,:,4)*dUx_dxi + Elem%Acoeff(:,:,:,12)*dUx_deta + Elem%Acoeff(:,:,:,19)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,25)*dUy_dxi + Elem%Acoeff(:,:,:,30)*dUy_deta + Elem%Acoeff(:,:,:,31)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,32)*dUz_dxi + Elem%Acoeff(:,:,:,33)*dUz_deta + Elem%Acoeff(:,:,:,34)*dUz_dzeta 

    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1.,t1(0,0,n_z), m1, htprimey ,m2, 0., s0(0,0,n_z),m1 )
    enddo
    Uyloc = s0 + Uyloc

    t1 = Elem%Acoeff(:,:,:,5)*dUx_dxi + Elem%Acoeff(:,:,:,13)*dUx_deta + Elem%Acoeff(:,:,:,20)* dUx_dzeta +&
         Elem%Acoeff(:,:,:,26)*dUy_dxi + Elem%Acoeff(:,:,:,31)*dUy_deta + Elem%Acoeff(:,:,:,35)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,36)*dUz_dxi + Elem%Acoeff(:,:,:,37)*dUz_deta + Elem%Acoeff(:,:,:,38)*dUz_dzeta

    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., t1(0,0,0), m1*m2, htprimez ,m3, 0., s0, m1*m2 )
    Uyloc = s0 + Uyloc

    t1 = Elem%Acoeff(:,:,:,6)*dUx_dxi + Elem%Acoeff(:,:,:,14)*dUx_deta + Elem%Acoeff(:,:,:,21)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,27)*dUy_dxi + Elem%Acoeff(:,:,:,32)*dUy_deta + Elem%Acoeff(:,:,:,36)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,39)*dUz_dxi + Elem%Acoeff(:,:,:,40)*dUz_deta + Elem%Acoeff(:,:,:,41)*dUz_dzeta

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,t1(0,0,0) ,m1, 0., Uzloc, m1 )

    t1 = Elem%Acoeff(:,:,:,7)*dUx_dxi + Elem%Acoeff(:,:,:,15)*dUx_deta + Elem%Acoeff(:,:,:,22)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,28)*dUy_dxi + Elem%Acoeff(:,:,:,33)*dUy_deta + Elem%Acoeff(:,:,:,37)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,40)*dUz_dxi + Elem%Acoeff(:,:,:,42)*dUz_deta + Elem%Acoeff(:,:,:,43)*dUz_dzeta

    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1.,t1(0,0,n_z), m1, htprimey ,m2, 0., s0(0,0,n_z),m1 )
    enddo
    Uzloc = s0 + Uzloc

    t1 = Elem%Acoeff(:,:,:,8)*dUx_dxi + Elem%Acoeff(:,:,:,16)*dUx_deta + Elem%Acoeff(:,:,:,23)*dUx_dzeta + &
         Elem%Acoeff(:,:,:,29)*dUy_dxi + Elem%Acoeff(:,:,:,34)*dUy_deta + Elem%Acoeff(:,:,:,38)*dUy_dzeta + &
         Elem%Acoeff(:,:,:,41)*dUz_dxi + Elem%Acoeff(:,:,:,43)*dUz_deta + Elem%Acoeff(:,:,:,44)*dUz_dzeta        

    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., t1(0,0,0), m1*m2, htprimez ,m3, 0., s0, m1*m2 )
    Uzloc = Uzloc + s0

    Elem%Forces(:,:,:,0) = Uxloc
    Elem%Forces(:,:,:,1) = Uyloc
    Elem%Forces(:,:,:,2) = Uzloc

    return
  end subroutine compute_InternalForces_Elem



  ! ###########################################################
  subroutine Prediction_Elem_PML_Veloc (Elem, bega, dt, hTprimex, Hprimey, Hprimez)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: bega, dt
    real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
    real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
    real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez

    real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
    real, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
    real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Vxloc,Vyloc,Vzloc, dVx_dxi,dVx_deta,dVx_dzeta, & 
         dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta

    integer :: m1, m2,m3 ,n_z,n_x,i,j,k


    m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

    Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2)  = Elem%Veloc(1:m1-2,1:m2-2,1:m3-2,0:2) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,0) ,m1, 0., dVx_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,0), m1, hprimey ,m2, 0., dVx_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0., dVx_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,1) ,m1, 0., dVy_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,1), m1, hprimey ,m2, 0., dVy_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0., dVy_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,2) ,m1, 0., dVz_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,2), m1, hprimey ,m2, 0., dVz_deta(0,0,n_z),m1)
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0., dVz_dzeta, m1*m2 )

    if (Elem%aniso) then

       Elem%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta &
            + Elem%Acoeff(:,:,:,2) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta &
            + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta &
            + Elem%Acoeff(:,:,:,8) * dVz_dzeta)
       Elem%Diagonal_Stress4a (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress4a (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta &
            + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
       Elem%Diagonal_Stress4b (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress4b (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta &
            + Elem%Acoeff(:,:,:,14) * dVy_dzeta)
       Elem%Diagonal_Stress5a (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress5a (:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVx_dxi + Elem%Acoeff(:,:,:,16) * dVx_deta &
            + Elem%Acoeff(:,:,:,17) * dVx_dzeta)
       Elem%Diagonal_Stress5b (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress5b (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta &
            + Elem%Acoeff(:,:,:,20) * dVz_dzeta)
       Elem%Diagonal_Stress6a (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress6a(:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVy_dxi + Elem%Acoeff(:,:,:,22) * dVy_deta &
            + Elem%Acoeff(:,:,:,23) * dVy_dzeta)
       Elem%Diagonal_Stress6b (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress6b(:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVz_dxi + Elem%Acoeff(:,:,:,25) * dVz_deta &
            + Elem%Acoeff(:,:,:,26) * dVz_dzeta)
       !---------------------------------------------------------------------
       Elem%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,27) * dVx_dxi + Elem%Acoeff(:,:,:,28) * dVx_deta &
            + Elem%Acoeff(:,:,:,29) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,30) * dVy_dxi + Elem%Acoeff(:,:,:,31) * dVy_deta &
            + Elem%Acoeff(:,:,:,32) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,33) * dVz_dxi + Elem%Acoeff(:,:,:,34) * dVz_deta &
            + Elem%Acoeff(:,:,:,35) * dVz_dzeta)
       Elem%Diagonal_Stress4a (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress4a (:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,36) * dVx_dxi + Elem%Acoeff(:,:,:,37) * dVx_deta &
            + Elem%Acoeff(:,:,:,38) * dVx_dzeta)
       Elem%Diagonal_Stress4b (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress4b (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,39) * dVy_dxi + Elem%Acoeff(:,:,:,40) * dVy_deta &
            + Elem%Acoeff(:,:,:,41) * dVy_dzeta)
       Elem%Diagonal_Stress5a (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress5a (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,42) * dVx_dxi + Elem%Acoeff(:,:,:,43) * dVx_deta &
            + Elem%Acoeff(:,:,:,44) * dVx_dzeta)
       Elem%Diagonal_Stress5b (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress5b(:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,45) * dVz_dxi + Elem%Acoeff(:,:,:,46) * dVz_deta &
            + Elem%Acoeff(:,:,:,47) * dVz_dzeta)
       Elem%Diagonal_Stress6a (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress6a (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,48) * dVy_dxi + Elem%Acoeff(:,:,:,49) * dVy_deta &
            + Elem%Acoeff(:,:,:,50) * dVy_dzeta)
       Elem%Diagonal_Stress6b (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress6b(:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,51) * dVz_dxi + Elem%Acoeff(:,:,:,52) * dVz_deta &
            + Elem%Acoeff(:,:,:,53) * dVz_dzeta)
       !---------------------------------------------------------------------
       Elem%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,54) * dVx_dxi + Elem%Acoeff(:,:,:,55) * dVx_deta &
            + Elem%Acoeff(:,:,:,56) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,57) * dVy_dxi + Elem%Acoeff(:,:,:,58) * dVy_deta &
            + Elem%Acoeff(:,:,:,59) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,60) * dVz_dxi + Elem%Acoeff(:,:,:,61) * dVz_deta &
            + Elem%Acoeff(:,:,:,62) * dVz_dzeta)
       Elem%Diagonal_Stress4a (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress4a (:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,63) * dVx_dxi + Elem%Acoeff(:,:,:,64) * dVx_deta &
            + Elem%Acoeff(:,:,:,65) * dVx_dzeta)
       Elem%Diagonal_Stress4b (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress4b(:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,66) * dVy_dxi + Elem%Acoeff(:,:,:,67) * dVy_deta &
            + Elem%Acoeff(:,:,:,68) * dVy_dzeta)
       Elem%Diagonal_Stress5a (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress5a(:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,69) * dVx_dxi + Elem%Acoeff(:,:,:,70) * dVx_deta &
            + Elem%Acoeff(:,:,:,71) * dVx_dzeta)
       Elem%Diagonal_Stress5b (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress5b(:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,72) * dVz_dxi + Elem%Acoeff(:,:,:,73) * dVz_deta &
            + Elem%Acoeff(:,:,:,74) * dVz_dzeta)
       Elem%Diagonal_Stress6a (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress6a(:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,75) * dVy_dxi + Elem%Acoeff(:,:,:,76) * dVy_deta &
            + Elem%Acoeff(:,:,:,77) * dVy_dzeta)
       Elem%Diagonal_Stress6b (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress6b(:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,78) * dVz_dxi + Elem%Acoeff(:,:,:,79) * dVz_deta &
            + Elem%Acoeff(:,:,:,80) * dVz_dzeta)
       !**********************************************************************
       Elem%Diagonal_Stress = Elem%Diagonal_Stress1 + Elem%Diagonal_Stress2 + Elem%Diagonal_Stress3 &
       + Elem%Diagonal_Stress4a + Elem%Diagonal_Stress4b + &
            Elem%Diagonal_Stress5a + Elem%Diagonal_Stress5b + Elem%Diagonal_Stress6a + Elem%Diagonal_Stress6b
       !**********************************************************************

       Elem%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,81) * dVx_dxi + Elem%Acoeff(:,:,:,82) * dVx_deta &
            + Elem%Acoeff(:,:,:,83) * dVx_dzeta)
       Elem%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,84) * dVy_dxi + Elem%Acoeff(:,:,:,85) * dVy_deta &
            + Elem%Acoeff(:,:,:,86) * dVy_dzeta)
       Elem%Residual_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress3 (:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,87) * dVz_dxi + Elem%Acoeff(:,:,:,88) * dVz_deta &
            + Elem%Acoeff(:,:,:,89) * dVz_dzeta)
       Elem%Residual_Stress4a (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress4a (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,90) * dVx_dxi + Elem%Acoeff(:,:,:,91) * dVx_deta &
            + Elem%Acoeff(:,:,:,92) * dVx_dzeta)
       Elem%Residual_Stress4b (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress4b (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,93) * dVy_dxi + Elem%Acoeff(:,:,:,94) * dVy_deta &
            + Elem%Acoeff(:,:,:,95) * dVy_dzeta)
       Elem%Residual_Stress5a (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress5a(:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,96) * dVx_dxi + Elem%Acoeff(:,:,:,97) * dVx_deta &
            + Elem%Acoeff(:,:,:,98) * dVx_dzeta)
       Elem%Residual_Stress5b (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress5b(:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,99) * dVz_dxi + Elem%Acoeff(:,:,:,100) * dVz_deta &
            + Elem%Acoeff(:,:,:,101) * dVz_dzeta)
       Elem%Residual_Stress6a (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress6a(:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,102) * dVy_dxi + Elem%Acoeff(:,:,:,103) * dVy_deta &
            + Elem%Acoeff(:,:,:,104) * dVy_dzeta)
       Elem%Residual_Stress6b (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress6b(:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,105) * dVz_dxi + Elem%Acoeff(:,:,:,106) * dVz_deta &
            + Elem%Acoeff(:,:,:,107) * dVz_dzeta)
       !---------------------------------------------------------------------
       Elem%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,108) * dVx_dxi + Elem%Acoeff(:,:,:,109) * dVx_deta &
            + Elem%Acoeff(:,:,:,110) * dVx_dzeta)
       Elem%Residual_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,111) * dVy_dxi + Elem%Acoeff(:,:,:,112) * dVy_deta &
            + Elem%Acoeff(:,:,:,113) * dVy_dzeta)
       Elem%Residual_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress3 (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,114) * dVz_dxi + Elem%Acoeff(:,:,:,115) * dVz_deta &
            + Elem%Acoeff(:,:,:,116) * dVz_dzeta)
       Elem%Residual_Stress4a (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress4a (:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,117) * dVx_dxi + Elem%Acoeff(:,:,:,118) * dVx_deta &
            + Elem%Acoeff(:,:,:,119) * dVx_dzeta)
       Elem%Residual_Stress4b (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress4b (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,120) * dVy_dxi + Elem%Acoeff(:,:,:,121) * dVy_deta &
            + Elem%Acoeff(:,:,:,122) * dVy_dzeta)
       Elem%Residual_Stress5a (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress5a(:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,123) * dVx_dxi + Elem%Acoeff(:,:,:,124) * dVx_deta &
            + Elem%Acoeff(:,:,:,125) * dVx_dzeta)
       Elem%Residual_Stress5b (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress5b (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,126) * dVz_dxi + Elem%Acoeff(:,:,:,127) * dVz_deta &
            + Elem%Acoeff(:,:,:,128) * dVz_dzeta)
       Elem%Residual_Stress6a (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress6a(:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,129) * dVy_dxi + Elem%Acoeff(:,:,:,130) * dVy_deta &
            + Elem%Acoeff(:,:,:,131) * dVy_dzeta)
       Elem%Residual_Stress6b (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress6b(:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,132) * dVz_dxi + Elem%Acoeff(:,:,:,133) * dVz_deta &
            + Elem%Acoeff(:,:,:,134) * dVz_dzeta)
       !---------------------------------------------------------------------
       Elem%Residual_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,135) * dVx_dxi + Elem%Acoeff(:,:,:,136) * dVx_deta &
            + Elem%Acoeff(:,:,:,137) * dVx_dzeta)
       Elem%Residual_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,138) * dVy_dxi + Elem%Acoeff(:,:,:,139) * dVy_deta &
            + Elem%Acoeff(:,:,:,140) * dVy_dzeta)
       Elem%Residual_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress3 (:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,141) * dVz_dxi + Elem%Acoeff(:,:,:,142) * dVz_deta &
            + Elem%Acoeff(:,:,:,143) * dVz_dzeta)
       Elem%Residual_Stress4a (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress4a(:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,144) * dVx_dxi + Elem%Acoeff(:,:,:,145) * dVx_deta &
            + Elem%Acoeff(:,:,:,146) * dVx_dzeta)
       Elem%Residual_Stress4b (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress4b(:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,147) * dVy_dxi + Elem%Acoeff(:,:,:,148) * dVy_deta &
            + Elem%Acoeff(:,:,:,149) * dVy_dzeta)
       Elem%Residual_Stress5a (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress5a(:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,150) * dVx_dxi + Elem%Acoeff(:,:,:,151) * dVx_deta &
            + Elem%Acoeff(:,:,:,152) * dVx_dzeta)
       Elem%Residual_Stress5b (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress5b(:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,153) * dVz_dxi + Elem%Acoeff(:,:,:,154) * dVz_deta &
            + Elem%Acoeff(:,:,:,155) * dVz_dzeta)
       Elem%Residual_Stress6a (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress6a(:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,156) * dVy_dxi + Elem%Acoeff(:,:,:,157) * dVy_deta &
            + Elem%Acoeff(:,:,:,158) * dVy_dzeta)
       Elem%Residual_Stress6b (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress6b(:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,159) * dVz_dxi + Elem%Acoeff(:,:,:,160) * dVz_deta &
            + Elem%Acoeff(:,:,:,161) * dVz_dzeta)
       !**********************************************************************
       Elem%Residual_Stress = Elem%Residual_Stress1 + Elem%Residual_Stress2 + Elem%Residual_Stress3 &
       + Elem%Residual_Stress4a + Elem%Residual_Stress4b &
       + Elem%Residual_Stress5a + Elem%Residual_Stress5b + Elem%Residual_Stress6a + Elem%Residual_Stress6b
       !**********************************************************************
    else 

       Elem%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta &
            + Elem%Acoeff(:,:,:,2) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta &
            + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,0) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta &
            + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

       Elem%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta &
            + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,1) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta &
            + Elem%Acoeff(:,:,:,14) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta &
            + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

       Elem%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,2) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta &
            + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
       Elem%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta &
            + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
       Elem%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta &
            + Elem%Acoeff(:,:,:,17) * dVz_dzeta)

       Elem%Diagonal_Stress = Elem%Diagonal_Stress1 + Elem%Diagonal_Stress2 + Elem%Diagonal_Stress3

       Elem%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,0) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta &
            + Elem%Acoeff(:,:,:,20) * dVy_dzeta)
       Elem%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,0) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta &
            + Elem%Acoeff(:,:,:,23) * dVx_dzeta)

       Elem%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,1) + &
            Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta &
            + Elem%Acoeff(:,:,:,20) * dVz_dzeta)
       Elem%Residual_Stress2 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,1) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta &
            + Elem%Acoeff(:,:,:,26) * dVx_dzeta)

       Elem%Residual_Stress1 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,2) + &
            Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta &
            + Elem%Acoeff(:,:,:,23) * dVz_dzeta)
       Elem%Residual_Stress2 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,2) + &
            Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta &
            + Elem%Acoeff(:,:,:,26) * dVy_dzeta)

       Elem%Residual_Stress = Elem%Residual_Stress1 + Elem%Residual_Stress2
    endif


    return
  end subroutine Prediction_Elem_PML_Veloc



  ! ###########################################################
  subroutine Prediction_Elem_FPML_Veloc (Elem, bega, dt, hTprimex, Hprimey, Hprimez,fil)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: bega, dt, fil
    real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
    real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
    real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez


    real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
    real, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
    real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Vxloc,Vyloc,Vzloc, dVx_dxi,dVx_deta,dVx_dzeta, & 
         dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta, Stress_ausiliar

    integer :: m1, m2,m3 ,n_z,n_x,i,j,k
    real :: fil2


    m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
    fil2 = fil**2

    Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2)  = Elem%Veloc(1:m1-2,1:m2-2,1:m3-2,0:2) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)


    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,0) ,m1, 0., dVx_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,0), m1, hprimey ,m2, 0., dVx_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0., dVx_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,1) ,m1, 0., dVy_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,1), m1, hprimey ,m2, 0., dVy_deta(0,0,n_z),m1 )
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0., dVy_dzeta, m1*m2 )

    call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,2) ,m1, 0., dVz_dxi, m1 )
    do n_z = 0,Elem%ngllz-1
       call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,2), m1, hprimey ,m2, 0., dVz_deta(0,0,n_z),m1)
    enddo
    call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0., dVz_dzeta, m1*m2 )


    Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,0)
    Elem%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,0) + &
         Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + &
         Elem%Acoeff(:,:,:,2) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,0)
    Elem%I_Diagonal_Stress1  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,0) )                  

    Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,0)
    Elem%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,0) + &
         Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + & 
         Elem%Acoeff(:,:,:,5) * dVy_dzeta) + Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,0)
    Elem%I_Diagonal_Stress2  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,0) )                  

    Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,0)
    Elem%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,0) + &
         Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + &
         Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,0)
    Elem%I_Diagonal_Stress3  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,0) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,0) )    

    Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,1)
    Elem%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,1) + &
         Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta &
         + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,1)
    Elem%I_Diagonal_Stress1  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,1) )                  

    Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,1)
    Elem%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,1) + &
         Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta &
         + Elem%Acoeff(:,:,:,14) * dVy_dzeta) + Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,1)
    Elem%I_Diagonal_Stress2  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,1) )                  

    Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,1)
    Elem%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,1) + &
         Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta &
         + Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,1)
    Elem%I_Diagonal_Stress3  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,1) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,1) )   

    Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,2)
    Elem%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,2) + &
         Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta &
         + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,2)
    Elem%I_Diagonal_Stress1  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,2) )                  

    Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,2)
    Elem%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,2) + &
         Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta &
         + Elem%Acoeff(:,:,:,5) * dVy_dzeta)+ Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,2)
    Elem%I_Diagonal_Stress2  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,2) )

    Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,2)
    Elem%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,2) + &
         Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta &
         + Elem%Acoeff(:,:,:,17) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,2)
    Elem%I_Diagonal_Stress3  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,2) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,2) )   


    Elem%Diagonal_Stress = Elem%Diagonal_Stress1 + Elem%Diagonal_Stress2 + Elem%Diagonal_Stress3


    Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,0)
    Elem%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,0) + &
         Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta &
         + Elem%Acoeff(:,:,:,20) * dVy_dzeta) + Elem%Isx * Elem%I_Residual_stress1(:,:,:,0)
    Elem%I_Residual_Stress1  (:,:,:,0)= Fil2* Elem%I_Residual_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,0) )             

    Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,0)
    Elem%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,0) + &
         Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta &
         + Elem%Acoeff(:,:,:,23) * dVx_dzeta) + Elem%Isy * Elem%I_Residual_stress2(:,:,:,0)
    Elem%I_Residual_Stress2  (:,:,:,0) = Fil2* Elem%I_Residual_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,0) ) 

    Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,1)
    Elem%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,1) + &
         Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta &
         + Elem%Acoeff(:,:,:,20) * dVz_dzeta) + Elem%Isx * Elem%I_Residual_stress1(:,:,:,1)
    Elem%I_Residual_Stress1  (:,:,:,1)= Fil2* Elem%I_Residual_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,1) )   

    Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,1)
    Elem%Residual_Stress2 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,1) + &
         Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta &
         + Elem%Acoeff(:,:,:,26) * dVx_dzeta) + Elem%Isz * Elem%I_Residual_stress2(:,:,:,1)
    Elem%I_Residual_Stress2  (:,:,:,1) = Fil2* Elem%I_Residual_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,1) ) 

    Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,2)
    Elem%Residual_Stress1 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,2) + &
         Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta &
         + Elem%Acoeff(:,:,:,23) * dVz_dzeta) + Elem%Isy * Elem%I_Residual_stress1(:,:,:,2)
    Elem%I_Residual_Stress1  (:,:,:,2)= Fil2* Elem%I_Residual_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,2) ) 

    Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,2)
    Elem%Residual_Stress2 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,2) + &
         Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta &
         + Elem%Acoeff(:,:,:,26) * dVy_dzeta)+ Elem%Isz * Elem%I_Residual_stress2(:,:,:,2)
    Elem%I_Residual_Stress2  (:,:,:,2) = Fil2* Elem%I_Residual_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * & 
         (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,2) ) 

    Elem%Residual_Stress = Elem%Residual_Stress1 + Elem%Residual_Stress2 

    return
  end subroutine Prediction_Elem_FPML_Veloc




  ! ###########################################################
  subroutine Correction_Elem_PML_Veloc (Elem, dt)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: dt

    integer :: ngllx, nglly, ngllz, i 


    ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz

    do i = 0,2
       Elem%Veloc1(:,:,:,i) = Elem%DumpVx(:,:,:,0) * Elem%Veloc1(:,:,:,i) + dt * &
       Elem%DumpVx(:,:,:,1)*Elem%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
       Elem%Veloc2(:,:,:,i) = Elem%DumpVy(:,:,:,0) * Elem%Veloc2(:,:,:,i) + dt * &
       Elem%DumpVy(:,:,:,1)*Elem%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
       Elem%Veloc3(:,:,:,i) = Elem%DumpVz(:,:,:,0) * Elem%Veloc3(:,:,:,i) + dt * &
       Elem%DumpVz(:,:,:,1)*Elem%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
    enddo

    Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Veloc1 + Elem%Veloc2 + Elem%Veloc3

    return
  end subroutine Correction_Elem_PML_Veloc

  ! ###########################################################
  subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, intent (IN) :: dt, fil

    integer :: ngllx, nglly, ngllz, i 
    real :: fil2
    real, dimension (1:Elem%ngllx-2,1:Elem%nglly-2,1:Elem%ngllz-2) :: Ausiliar_velocity

    ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz
    fil2 = fil**2

    do i = 0,2
       Ausiliar_Velocity =  Elem%Veloc1(:,:,:,i)
       Elem%Veloc1(:,:,:,i) = Elem%DumpVx(:,:,:,0) * Elem%Veloc1(:,:,:,i) + dt * &
       Elem%DumpVx(:,:,:,1)*Elem%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivx * Elem%Iveloc1 (:,:,:,i)
       Elem%Iveloc1 (:,:,:,i)  = Fil2*Elem%Iveloc1 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc1(:,:,:,i))

       Ausiliar_Velocity =  Elem%Veloc2(:,:,:,i)
       Elem%Veloc2(:,:,:,i) = Elem%DumpVy(:,:,:,0) * Elem%Veloc2(:,:,:,i) + dt * &
       Elem%DumpVy(:,:,:,1)*Elem%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivy * Elem%Iveloc2 (:,:,:,i)
       Elem%Iveloc2 (:,:,:,i)  = Fil2*Elem%Iveloc2 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc2(:,:,:,i))

       Ausiliar_Velocity =  Elem%Veloc3(:,:,:,i)
       Elem%Veloc3(:,:,:,i) = Elem%DumpVz(:,:,:,0) * Elem%Veloc3(:,:,:,i) + dt * &
       Elem%DumpVz(:,:,:,1)*Elem%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivz * Elem%Iveloc3 (:,:,:,i)
       Elem%Iveloc3 (:,:,:,i)  = Fil2*Elem%Iveloc3 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc3(:,:,:,i))
    enddo

    Elem%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Veloc1 + Elem%Veloc2 + Elem%Veloc3

    return
  end subroutine Correction_Elem_FPML_Veloc


!!$! ###########################################################
!!$subroutine  compute_InternalForces_PML_Elem (Elem, hprimex, hTprimey,htprimez)
!!$
!!$  implicit none
!!$
!!$  type (Element), intent (INOUT) :: Elem
!!$  real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
!!$  real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
!!$  real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
!!$
!!$  integer :: m1, m2, m3, n_z
!!$  real, dimension ( 0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1)  :: s0,s1
!!$
!!$
!!$  m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
!!$
!!$  s0 = Elem%Acoeff(:,:,:,27) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,0) + & 
!!$       Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,1)
!!$  call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
!!$  Elem%Forces1(:,:,:,0) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,30) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,0) + & 
!!$       Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,1)
!!$  do n_z = 0,m3-1
!!$     call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
!!$  enddo
!!$  Elem%Forces2(:,:,:,0) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,33) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,0) + & 
!!$       Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,1)
!!$
!!$  call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
!!$  Elem%Forces3(:,:,:,0) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Diagonal_Stress(:,:,:,1) + & 
!!$       Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,2)
!!$  call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
!!$  Elem%Forces1(:,:,:,1) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Diagonal_Stress(:,:,:,1) + & 
!!$       Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,2)
!!$  do n_z = 0,m3-1
!!$     call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
!!$  enddo
!!$  Elem%Forces2(:,:,:,1) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Diagonal_Stress(:,:,:,1) + &
!!$       Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,2)
!!$
!!$  call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
!!$  Elem%Forces3(:,:,:,1) = s1
!!$
!!$
!!$  s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,2) + & 
!!$       Elem%Acoeff(:,:,:,29) * Elem%Diagonal_Stress(:,:,:,2)
!!$  call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
!!$  Elem%Forces1(:,:,:,2) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,2) + & 
!!$       Elem%Acoeff(:,:,:,32) * Elem%Diagonal_Stress(:,:,:,2)
!!$  do n_z = 0,m3-1
!!$     call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
!!$  enddo
!!$  Elem%Forces2(:,:,:,2) = s1
!!$
!!$  s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,2) + & 
!!$       Elem%Acoeff(:,:,:,35) * Elem%Diagonal_Stress(:,:,:,2)
!!$
!!$  call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
!!$  Elem%Forces3(:,:,:,2) = s1
!!$
!!$  Elem%Forces = Elem%Forces1 + Elem%Forces2 + Elem%Forces3
!!$
!!$  return
!!$end subroutine compute_InternalForces_PML_Elem




  ! ###########################################################
  subroutine  compute_InternalForces_PML_Elem (Elem, hprimex, hTprimey,htprimez)

    implicit none

    type (Element), intent (INOUT) :: Elem
    real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
    real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
    real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez

    integer :: m1, m2, m3, n_z
    real, dimension ( 0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1)  :: s0,s1


    m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

    if (Elem%aniso) then
       s0 = Elem%Acoeff(:,:,:,162) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,163) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,164) * Elem%Residual_Stress(:,:,:,1)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,165) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,166) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,167) * Elem%Residual_Stress(:,:,:,1)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,168) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,169) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,170) * Elem%Residual_Stress(:,:,:,1)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,162) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,163) * Elem%Diagonal_Stress(:,:,:,1) + & 
            Elem%Acoeff(:,:,:,164) * Elem%Residual_Stress(:,:,:,2)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,1) = s1

       s0 = Elem%Acoeff(:,:,:,165) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,166) * Elem%Diagonal_Stress(:,:,:,1) + & 
            Elem%Acoeff(:,:,:,167) * Elem%Residual_Stress(:,:,:,2)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,1) = s1

       s0 = Elem%Acoeff(:,:,:,168) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,169) * Elem%Diagonal_Stress(:,:,:,1) + &
            Elem%Acoeff(:,:,:,170) * Elem%Residual_Stress(:,:,:,2)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,1) = s1


       s0 = Elem%Acoeff(:,:,:,162) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,163) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,164) * Elem%Diagonal_Stress(:,:,:,2)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,2) = s1

       s0 = Elem%Acoeff(:,:,:,165) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,166) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,167) * Elem%Diagonal_Stress(:,:,:,2)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,2) = s1

       s0 = Elem%Acoeff(:,:,:,168) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,169) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,170) * Elem%Diagonal_Stress(:,:,:,2)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,2) = s1

    else
       s0 = Elem%Acoeff(:,:,:,27) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,1)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,30) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,1)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,33) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,0) + & 
            Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,1)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,0) = s1

       s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Diagonal_Stress(:,:,:,1) + & 
            Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,2)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,1) = s1

       s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Diagonal_Stress(:,:,:,1) + & 
            Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,2)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,1) = s1

       s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Diagonal_Stress(:,:,:,1) + &
            Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,2)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,1) = s1


       s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,29) * Elem%Diagonal_Stress(:,:,:,2)
       call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
       Elem%Forces1(:,:,:,2) = s1

       s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,32) * Elem%Diagonal_Stress(:,:,:,2)
       do n_z = 0,m3-1
          call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
       enddo
       Elem%Forces2(:,:,:,2) = s1

       s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,2) + & 
            Elem%Acoeff(:,:,:,35) * Elem%Diagonal_Stress(:,:,:,2)

       call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
       Elem%Forces3(:,:,:,2) = s1
    endif


    Elem%Forces = Elem%Forces1 + Elem%Forces2 + Elem%Forces3

    return
  end subroutine compute_InternalForces_PML_Elem

  ! ###########################################################
end module selement

