subroutine allocate_domain (Tdomain,rg)

    ! Modified by Gaetano Festa 23/2/2005
    ! Modified by Paul Cupillard 09/12/2005


    use sdomain

    implicit none

    type(domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: rg

    integer :: n, nf,ne,nv, i,j,k, ngllx,nglly,ngllz, ngll1,ngll2, ngll, ngllPML, ngllSO, mat,ns,n_point_face
    real::x0,x1,y0,y3,z0,z5
    logical::YesInSlice

    !if (Tdomain%Field_Order(6)) Tdomain%MaxiSliceSize=50000
    do n=0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        Tdomain%specel(n)%aniso = .false.
        if (Tdomain%sSubDomain(mat)%material_type(2:2) == "A") Tdomain%specel(n)%aniso = .true.



        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        allocate (Tdomain%specel(n)%Density (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        if (.not. Tdomain%specel(n)%aniso) then
            allocate (Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
            allocate (Tdomain%specel(n)%Mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        endif
        allocate (Tdomain%specel(n)%MassMat (0:ngllx-1, 0:nglly-1, 0:ngllz-1))

        if (Tdomain%TimeD%velocity_scheme) then
            !!$        allocate (Tdomain%specel(n)%Veloc (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%Veloc (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%Accel (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%V0 (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            allocate (Tdomain%specel(n)%Forces( 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
            !!$if (Tdomain%TimeReverse) allocate (Tdomain%specel(n)%D0 (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
            Tdomain%specel(n)%Veloc = 0
            Tdomain%specel(n)%Accel = 0
            Tdomain%specel(n)%V0 = 0
            Tdomain%specel(n)%Forces = 0
            if (Tdomain%specel(n)%PML) then
                if (Tdomain%specel(n)%aniso) then
                    allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:170))
                    allocate (Tdomain%specel(n)%Diagonal_Stress4a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Diagonal_Stress4b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Diagonal_Stress5a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Diagonal_Stress5b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Diagonal_Stress6a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Diagonal_Stress6b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))

                    allocate (Tdomain%specel(n)%Residual_Stress3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress4a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress4b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress5a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress5b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress6a (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                    allocate (Tdomain%specel(n)%Residual_Stress6b (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))

                    Tdomain%specel(n)%Diagonal_Stress4a = 0.
                    Tdomain%specel(n)%Diagonal_Stress4b = 0.
                    Tdomain%specel(n)%Diagonal_Stress5a = 0.
                    Tdomain%specel(n)%Diagonal_Stress5b = 0.
                    Tdomain%specel(n)%Diagonal_Stress6a = 0.
                    Tdomain%specel(n)%Diagonal_Stress6b = 0.

                    Tdomain%specel(n)%Residual_Stress3 = 0.
                    Tdomain%specel(n)%Residual_Stress4a = 0.
                    Tdomain%specel(n)%Residual_Stress4b = 0.
                    Tdomain%specel(n)%Residual_Stress5a = 0.
                    Tdomain%specel(n)%Residual_Stress5b = 0.
                    Tdomain%specel(n)%Residual_Stress6a = 0.
                    Tdomain%specel(n)%Residual_Stress6b = 0.
                else
                    allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:35))
                endif
                allocate (Tdomain%specel(n)%Diagonal_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Diagonal_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Diagonal_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Diagonal_Stress3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Residual_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Residual_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Residual_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Veloc1 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
                allocate (Tdomain%specel(n)%Veloc2 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
                allocate (Tdomain%specel(n)%Veloc3 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
                allocate (Tdomain%specel(n)%Forces1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Forces2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%Forces3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%DumpSx (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
                allocate (Tdomain%specel(n)%DumpSy (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
                allocate (Tdomain%specel(n)%DumpSz (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
                allocate (Tdomain%specel(n)%DumpMass (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                allocate (Tdomain%specel(n)%DumpVx (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
                allocate (Tdomain%specel(n)%DumpVy (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
                allocate (Tdomain%specel(n)%DumpVz (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
                Tdomain%specel(n)%Diagonal_Stress = 0
                Tdomain%specel(n)%Diagonal_Stress1 = 0
                Tdomain%specel(n)%Diagonal_Stress2 = 0
                Tdomain%specel(n)%Diagonal_Stress3 = 0
                Tdomain%specel(n)%Residual_Stress = 0
                Tdomain%specel(n)%Residual_Stress1 = 0
                Tdomain%specel(n)%Residual_Stress2 = 0
                Tdomain%specel(n)%Veloc1 = 0
                Tdomain%specel(n)%Veloc2 = 0
                Tdomain%specel(n)%Veloc3 = 0
                if (Tdomain%specel(n)%FPML) then
                    allocate (Tdomain%specel(n)%Isx(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%Isy(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%Isz(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%Ivx(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%Ivy(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%Ivz(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate (Tdomain%specel(n)%I_Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate (Tdomain%specel(n)%I_Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate (Tdomain%specel(n)%I_Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate (Tdomain%specel(n)%I_Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate (Tdomain%specel(n)%I_Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate (Tdomain%specel(n)%IVeloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    allocate (Tdomain%specel(n)%IVeloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    allocate (Tdomain%specel(n)%IVeloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    Tdomain%specel(n)%I_Diagonal_Stress1 = 0.
                    Tdomain%specel(n)%I_Diagonal_Stress2 = 0.
                    Tdomain%specel(n)%I_Diagonal_Stress3 = 0.
                    Tdomain%specel(n)%I_Residual_Stress1 = 0.
                    Tdomain%specel(n)%I_Residual_Stress2 = 0.
                    Tdomain%specel(n)%IVeloc1 = 0.
                    Tdomain%specel(n)%IVeloc2 = 0.
                    Tdomain%specel(n)%IVeloc3 = 0.
                endif
            else
                allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:44))
                !!$           allocate (Tdomain%specel(n)%Displ (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
                allocate (Tdomain%specel(n)%Displ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
                Tdomain%specel(n)%Displ = 0.
                if (((Tdomain%logicD%save_snapshots.eqv..true.) &
                .and.  (Tdomain%Field_Order(2).eqv..true.)) .or. (Tdomain%logicD%save_energy.eqv..true.)) then
                    allocate(Tdomain%specel(n)%c(1:6,1:6,0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                !!$              allocate(Tdomain%specel(n)%stress_sigma(1:6,0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                endif
                if ((Tdomain%logicD%save_snapshots.eqv..true.) &
                .and.  (Tdomain%Field_Order(4).eqv..true.)) then
                    allocate (Tdomain%specel(n)%MaxAbsVeloc(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
!!!                    allocate (Tdomain%specel(n)&
!!!!                    %TravelTime (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 1:size(Tdomain%Ratios)+1))
                    allocate (Tdomain%specel(n)%TravelTimeFound (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 1:size(Tdomain%Ratios)))
                    Tdomain%specel(n)%MaxAbsVeloc=0.
!!!                    Tdomain%specel(n)%TravelTime=0.
                    Tdomain%specel(n)%TravelTimeFound=.false.
                endif

            endif
        endif
        if (Tdomain%Field_Order(6) .or. Tdomain%Field_Order(7)) then
            x0=Tdomain%Coord_nodes(0,Tdomain%specel(n)%Control_Nodes(0))
            x1=Tdomain%Coord_nodes(0,Tdomain%specel(n)%Control_Nodes(1))

            y0=Tdomain%Coord_nodes(1,Tdomain%specel(n)%Control_Nodes(0))
            y3=Tdomain%Coord_nodes(1,Tdomain%specel(n)%Control_Nodes(3))
        
            z0=Tdomain%Coord_nodes(2,Tdomain%specel(n)%Control_Nodes(0))
            z5=Tdomain%Coord_nodes(2,Tdomain%specel(n)%Control_Nodes(5))


            do ns=1,Tdomain%nbSlices

                select case(int(Tdomain%Slices(ns,1)))
                    case (0)
                        YesInSlice=(abs(x0-Tdomain%Slices(ns,2)) .lt. 1e-4) .and. &
                            ((y0-Tdomain%Slices(ns,3)) .gt. -1e-4 ) .and. &
                            ((y3-Tdomain%Slices(ns,4)) .lt.  1e-4 ) .and. &
                            ((z0-Tdomain%Slices(ns,5)) .gt. -1e-4 ) .and. &
                            ((z5-Tdomain%Slices(ns,6)) .lt.  1e-4 )
                        n_point_face=nglly*ngllz
                    !if (abs(x0-Tdomain%Slices(ns,2)) .lt. 1e-4) print *, y0-Tdomain%Slices(ns,3),y3-Tdomain%Slices(ns,4),z0-Tdomain%Slices(ns,5),z5-Tdomain%Slices(ns,6),'rank',rg
                    case (1)
                        YesInSlice=(abs(y0-Tdomain%Slices(ns,2)) .lt. 1e-4) .and. &
                            ((x0-Tdomain%Slices(ns,3)) .gt. -1e-4 ) .and. &
                            ((x1-Tdomain%Slices(ns,4)) .lt.  1e-4 ) .and. &
                            ((z0-Tdomain%Slices(ns,5)) .gt. -1e-4 ) .and. &
                            ((z5-Tdomain%Slices(ns,6)) .lt.  1e-4 )
                        n_point_face=ngllx*ngllz
                    case (2)
                        YesInSlice=(abs(z5-Tdomain%Slices(ns,2)) .lt. 1e-4) .and. &
                            ((x0-Tdomain%Slices(ns,3)) .gt. -1e-4 ) .and. &
                            ((x1-Tdomain%Slices(ns,4)) .lt.  1e-4 ) .and. &
                            ((y0-Tdomain%Slices(ns,5)) .gt. -1e-4 ) .and. &
                            ((y3-Tdomain%Slices(ns,6)) .lt.  1e-4 )
                        n_point_face=ngllx*nglly
                endselect

                if (YesInSlice) then
                    Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize-1)=Tdomain&
                    %I_elem_inSlices(ns,Tdomain%MaxiSliceSize-1)+n_point_face
                    Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize)=Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize)+1
                    Tdomain%I_elem_inSlices(ns,Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize))=n
                endif

            enddo
        endif
    enddo
    !print *, Tdomain%I_elem_inSlices,size(Tdomain%I_elem_inSlices),rg
    do n = 0, Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        allocate (Tdomain%sFace(n)%MassMat (1:ngll1-2, 1:ngll2-2))
        allocate (Tdomain%sFace(n)%Veloc (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%Forces (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%Accel (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%V0 (1:ngll1-2, 1:ngll2-2, 0:2))
        !!$if (Tdomain%TimeReverse) allocate (Tdomain%sFace(n)%D0 (1:ngll1-2, 1:ngll2-2, 0:2))
        Tdomain%sFace(n)%MassMat = 0
        Tdomain%sFace(n)%Veloc = 0
        Tdomain%sFace(n)%Accel = 0
        Tdomain%sFace(n)%V0 = 0
        Tdomain%sFace(n)%Forces = 0
        if (Tdomain%sFace(n)%PML ) then
            allocate (Tdomain%sFace(n)%Forces1 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%Forces2 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%Forces3 (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%DumpMass (1:ngll1-2, 1:ngll2-2, 0:2))
            allocate (Tdomain%sFace(n)%Veloc1 (1:ngll1-2, 1:ngll2-2,0:2))
            allocate (Tdomain%sFace(n)%Veloc2 (1:ngll1-2, 1:ngll2-2,0:2))
            allocate (Tdomain%sFace(n)%Veloc3 (1:ngll1-2, 1:ngll2-2,0:2))
            allocate (Tdomain%sFace(n)%DumpVx (1:ngll1-2, 1:ngll2-2, 0:1))
            allocate (Tdomain%sFace(n)%DumpVy (1:ngll1-2, 1:ngll2-2, 0:1))
            allocate (Tdomain%sFace(n)%DumpVz (1:ngll1-2, 1:ngll2-2, 0:1))
            Tdomain%sFace(n)%DumpMass = 0.
            Tdomain%sFace(n)%Veloc1 = 0.
            Tdomain%sFace(n)%Veloc2 = 0.
            Tdomain%sFace(n)%Veloc3 = 0.
            if (Tdomain%sFace(n)%FPML ) then
                allocate (Tdomain%sFace(n)%Ivx(1:ngll1-2,1:ngll2-1))
                allocate (Tdomain%sFace(n)%Ivy(1:ngll1-2,1:ngll2-1))
                allocate (Tdomain%sFace(n)%Ivz(1:ngll1-2,1:ngll2-1))
                allocate (Tdomain%sFace(n)%IVeloc1(1:ngll1-2,1:ngll2-1,0:2))
                allocate (Tdomain%sFace(n)%IVeloc2(1:ngll1-2,1:ngll2-1,0:2))
                allocate (Tdomain%sFace(n)%IVeloc3(1:ngll1-2,1:ngll2-1,0:2))
                Tdomain%sFace(n)%IVeloc1 = 0.
                Tdomain%sFace(n)%IVeloc2 = 0.
                Tdomain%sFace(n)%IVeloc3 = 0.
                Tdomain%sFace(n)%Ivx = 0.
                Tdomain%sFace(n)%Ivy = 0.
                Tdomain%sFace(n)%Ivz = 0.
            endif
        else
            allocate (Tdomain%sFace(n)%Displ (1:ngll1-2, 1:ngll2-2, 0:2))
            Tdomain%sFace(n)%Displ = 0
        endif
    enddo

    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        allocate (Tdomain%sEdge(n)%MassMat (1:ngll-2))
        allocate (Tdomain%sEdge(n)%Veloc (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Forces (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Accel (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%V0 (1:ngll-2, 0:2))
        !!$ if (Tdomain%TimeReverse) allocate (Tdomain%sEdge(n)%D0 (1:ngll-2, 0:2))
        Tdomain%sEdge(n)%MassMat = 0
        Tdomain%sEdge(n)%Veloc = 0
        Tdomain%sEdge(n)%Accel = 0
        Tdomain%sEdge(n)%V0 = 0
        Tdomain%sEdge(n)%Forces = 0
        if (Tdomain%sEdge(n)%PML) then
            allocate (Tdomain%sEdge(n)%Forces1 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%Forces2 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%Forces3 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%DumpMass (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%Veloc1 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%Veloc2 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%Veloc3 (1:ngll-2, 0:2))
            allocate (Tdomain%sEdge(n)%DumpVx (1:ngll-2, 0:1))
            allocate (Tdomain%sEdge(n)%DumpVy (1:ngll-2, 0:1))
            allocate (Tdomain%sEdge(n)%DumpVz (1:ngll-2, 0:1))
            Tdomain%sEdge(n)%DumpMass = 0
            Tdomain%sEdge(n)%Veloc1 = 0
            Tdomain%sEdge(n)%Veloc2 = 0
            Tdomain%sEdge(n)%Veloc3 = 0
            if (Tdomain%sEdge(n)%FPML ) then
                allocate (Tdomain%sEdge(n)%Ivx(1:ngll-2))
                allocate (Tdomain%sEdge(n)%Ivy(1:ngll-2))
                allocate (Tdomain%sEdge(n)%Ivz(1:ngll-2))
                allocate (Tdomain%sEdge(n)%IVeloc1(1:ngll-2,0:2))
                allocate (Tdomain%sEdge(n)%IVeloc2(1:ngll-2,0:2))
                allocate (Tdomain%sEdge(n)%IVeloc3(1:ngll-2,0:2))
                Tdomain%sEdge(n)%IVeloc1 = 0.
                Tdomain%sEdge(n)%IVeloc2 = 0.
                Tdomain%sEdge(n)%IVeloc3 = 0.
                Tdomain%sEdge(n)%Ivx = 0.
                Tdomain%sEdge(n)%Ivy = 0.
                Tdomain%sEdge(n)%Ivz = 0.
            endif
        else
            allocate (Tdomain%sEdge(n)%Displ (1:ngll-2, 0:2))
            Tdomain%sEdge(n)%Displ = 0
        endif
    enddo

    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%MassMat = 0
        Tdomain%sVertex(n)%Veloc = 0
        Tdomain%sVertex(n)%Accel = 0
        Tdomain%sVertex(n)%V0 = 0
        Tdomain%sVertex(n)%Forces = 0
        !!$Tdomain%sVertex(n)%D0 = 0
        if (Tdomain%sVertex(n)%PML) then
            Tdomain%sVertex(n)%DumpMass = 0
            Tdomain%sVertex(n)%Veloc1 = 0
            Tdomain%sVertex(n)%Veloc2 = 0
            Tdomain%sVertex(n)%Veloc3 = 0
            if (Tdomain%sVertex(n)%FPML) then
                allocate (Tdomain%sVertex(n)%IVeloc1(0:2))
                allocate (Tdomain%sVertex(n)%IVeloc2(0:2))
                allocate (Tdomain%sVertex(n)%IVeloc3(0:2))
                allocate (Tdomain%sVertex(n)%Ivx(0:0))
                allocate (Tdomain%sVertex(n)%Ivy(0:0))
                allocate (Tdomain%sVertex(n)%Ivz(0:0))
                Tdomain%sVertex(n)%IVeloc1 = 0
                Tdomain%sVertex(n)%IVeloc2 = 0
                Tdomain%sVertex(n)%IVeloc3 = 0
                Tdomain%sVertex(n)%Ivx= 0
                Tdomain%sVertex(n)%Ivy= 0
                Tdomain%sVertex(n)%Ivz= 0
            endif
        else
            Tdomain%sVertex(n)%Displ = 0
        endif
    enddo

    if ( Tdomain%n_proc > 1 ) then
        do n = 0,Tdomain%n_proc-1

            ngll = 0
            ngllPML = 0
            ngllSO = 0
            do i = 0,Tdomain%sComm(n)%nb_faces-1
                nf = Tdomain%sComm(n)%faces(i)
                ngll1 = Tdomain%sFace(nf)%ngll1; ngll2 = Tdomain%sFace(nf)%ngll2
                ngll = ngll + (ngll1-2)*(ngll2-2)
                if (Tdomain%sFace(nf)%PML) then
                    ngllPML = ngllPML + (ngll1-2)*(ngll2-2)
                   !if (rg==25 .and. n==27) print*,'alloface',rg,n,ngllPML
                   !if (rg==27 .and. n==25) write(60,*)'ngllPML27face',ngllPML
                   !if (rg==25 .and. n==27) print*,'ngllPML27face',ngllPML
                endif
            enddo
            do i = 0,Tdomain%sComm(n)%nb_edges-1
                ne = Tdomain%sComm(n)%edges(i)
                ngll = ngll + Tdomain%sEdge(ne)%ngll-2
                if (Tdomain%sEdge(ne)%PML)  then
                    !if (rg==27 .and. n==25) write(60,*) ne
                    !if (rg==25 .and. n==27) print*, ne
                    !if (rg==27 .and. n==25) write(60,*)'ngllPML27edge',ngllPML
                    !if (rg==25 .and. n==27) print*,'ngllPML27edge',ngllPML
                    !if (rg==25 .and. n==27) print*,'alloedge',rg,n,ngllPML
                    ngllPML = ngllPML + Tdomain%sEdge(ne)%ngll-2
                   !if (rg==27 .and. n==25) write(60,*)'ngllPML27edge',ngllPML
                   !if (rg==25 .and. n==27) print*,'ngllPML27edge',ngllPML
                   !if (rg==25 .and. n==27) print*,'alloedge',rg,n,ngllPML
                endif
            enddo
            do i = 0,Tdomain%sComm(n)%nb_vertices-1
                nv = Tdomain%sComm(n)%vertices(i)
                ngll = ngll + 1
                if (Tdomain%sVertex(nv)%PML) then
                    !if (rg==27 .and. n==25) write(60,*) nv
                    !if (rg==25 .and. n==27) print*, nv
                    !if (rg==27 .and. n==25) write(60,*)'ngllPML27vert',ngllPML
                    !if (rg==25 .and. n==27) print*,'ngllPML27vert',ngllPML
                    !if (rg==25 .and. n==27) print*,'allovert',rg,n,ngllPML
                    ngllPML = ngllPML + 1
                   !if (rg==27 .and. n==25) write(60,*)'ngllPML27vert',ngllPML
                   !if (rg==25 .and. n==27) print*,'ngllPML27vert',ngllPML
                   !if (rg==25 .and. n==27) print*,'allovert',rg,n,ngllPML
                endif
            enddo
            do i = 0,Tdomain%sComm(n)%nb_edges_so-1
                ne = Tdomain%sComm(n)%edges_SO(i)
                ngllSO = ngllSO + Tdomain%sPlaneW%pEdge(ne)%ngll-2
            enddo
            do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
                ngllSO = ngllSO + 1
            enddo
            do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
                ne = Tdomain%sComm(n)%edges_Neu(i)
                ngllSO = ngllSO + Tdomain%sNeu%nEdge(ne)%ngll-2
            enddo
            do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
                ngllSO = ngllSO + 1
            enddo
            if (ngll>0) then
                allocate (Tdomain%sComm(n)%Give (0:ngll-1))
                allocate (Tdomain%sComm(n)%Take (0:ngll-1))
                allocate (Tdomain%sComm(n)%GiveForces (0:ngll-1, 0:2))
                allocate (Tdomain%sComm(n)%TakeForces (0:ngll-1, 0:2))
            endif
            if (ngllPML>0) then
                if (Tdomain%any_FPML ) then
                    allocate (Tdomain%sComm(n)%GivePML (0:ngllPML-1, 0:5))
                    allocate (Tdomain%sComm(n)%TakePML (0:ngllPML-1, 0:5))
                else
                    allocate (Tdomain%sComm(n)%GivePML (0:ngllPML-1, 0:2))
                    allocate (Tdomain%sComm(n)%TakePML (0:ngllPML-1, 0:2))
                endif
                allocate (Tdomain%sComm(n)%GiveForcesPML (0:ngllPML-1, 1:3, 0:2))
                allocate (Tdomain%sComm(n)%TakeForcesPML (0:ngllPML-1, 1:3, 0:2))
            endif
            if (ngllSO>0) then
                allocate (Tdomain%sComm(n)%GiveSO (0:ngllSO-1,0:2))
                allocate (Tdomain%sComm(n)%TakeSO (0:ngllSO-1,0:2))
            endif
            Tdomain%sComm(n)%ngll = ngll
            Tdomain%sComm(n)%ngllPML = ngllPML
            Tdomain%sComm(n)%ngllSO = ngllSO
           !if (rg==27) print*,'allo',rg,n,Tdomain%sComm(n)%ngll,Tdomain%sComm(n)%ngllPML,Tdomain%sComm(n)%ngllSO
           !if (n==27) print*,'allo',rg,n,Tdomain%sComm(n)%ngll,Tdomain%sComm(n)%ngllPML,Tdomain%sComm(n)%ngllSO
           !if (rg==25 .and. n==27) print*,'allo',rg,n,Tdomain%sComm(n)%ngll,Tdomain%sComm(n)%ngllPML,Tdomain%sComm(n)%ngllSO,ngllPML
           ! if (rg==27 .and. n==25) print*,'allo27-25',Tdomain%sComm(n)%nb_faces,Tdomain%sComm(n)%nb_edges,&
           !                                       Tdomain%sComm(n)%nb_vertices,Tdomain%sComm(n)%nb_edges_so,&
           !                                       Tdomain%sComm(n)%nb_vertices_so,Tdomain%sComm(n)%nb_edges_neu,&
           !                                       Tdomain%sComm(n)%nb_vertices_neu
           ! if (rg==25 .and. n==27) print*,'allo25-27',Tdomain%sComm(n)%nb_faces,Tdomain%sComm(n)%nb_edges,&
           !                                       Tdomain%sComm(n)%nb_vertices,Tdomain%sComm(n)%nb_edges_so,&
           !                                       Tdomain%sComm(n)%nb_vertices_so,Tdomain%sComm(n)%nb_edges_neu,&
           !                                       Tdomain%sComm(n)%nb_vertices_neu
        enddo
    endif

    !do n = 0, Tdomain%n_receivers-1
    !   allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%TimeD%NtimeMax-1, 0:2) )
    !enddo

    !print*,Tdomain%sComm(n)%nb_faces,Tdomain%sComm(n)%nb_edges,

    return
end subroutine allocate_domain
