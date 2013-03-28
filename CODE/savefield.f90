subroutine savefield (Tdomain,it,rg,icount)

  ! Modified by Elise Delavaud (26/02/06)


  use sdomain

  implicit none

  include 'mpif.h'


  type (domain), intent (INOUT):: Tdomain
  integer, intent (IN) :: it,rg,icount

  ! local variables
  integer :: i,n,nv,nbvert,kcount,nv_aus
  integer, dimension(:), allocatable:: count
  character (len=100) :: fnamef

  integer, dimension(:),allocatable ::NbElem, NbGlob 
  integer::SZG,SZG6,n_x,n_y,n_z,i1,NbOctInt, NbOctReal,desc,code,nbprocs,ngllx,nglly,ngllz,mat,l

  real, dimension(:,:,:), allocatable ::  dUx_dxi, dUx_deta, dUx_dzeta, &
       dUy_dxi, dUy_deta, dUy_dzeta, &
       dUz_dxi, dUz_deta, dUz_dzeta, &
       xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz
  integer,dimension(:,:,:),allocatable::GlobId
  real, dimension (:,:,:,:), allocatable :: stress_sigma,strain,Veloc_GLLs

  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer (kind=MPI_OFFSET_KIND) :: PosFile, PosFileI, PosFileR
  real :: AbsVeloc


  if (Tdomain%Field_Order(1)) then !Veloc Order only at vertex
     write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldx",icount
     open (61, file=fnamef, status="unknown", form="formatted")
     write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldy",icount
     open (62, file=fnamef, status="unknown", form="formatted")
     write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldz",icount
     open (63, file=fnamef, status="unknown", form="formatted")

     allocate (count(0:Tdomain%n_vertex-1))
     count = -1

     kcount = 0
     do n = 0, Tdomain%n_elem - 1
        !  if ( Tdomain%specel(n)%mat_index == 0 .or. Tdomain%specel(n)%mat_index == 1) then
        do i = 0,7
           nv = Tdomain%specel(n)%Near_Vertices(i)
           if ( count(nv) < 0 ) then
              kcount = kcount+1
              count(nv) = 1
           endif
        enddo
        !  endif
     enddo

     write (61,*) kcount,Tdomain%TimeD%time_snapshots_step
     write (62,*) kcount,Tdomain%TimeD%time_snapshots_step
     write (63,*) kcount,Tdomain%TimeD%time_snapshots_step

     count = -1
     do n = 0, Tdomain%n_elem - 1
        !  if ( Tdomain%specel(n)%mat_index == 0 .or. Tdomain%specel(n)%mat_index == 1 ) then
        do i = 0,7
           nv = Tdomain%specel(n)%Near_Vertices(i)
           if ( count(nv) < 0 ) then
              write (61,*) nv+1, Tdomain%svertex(nv)%Veloc(0)
              write (62,*) nv+1, Tdomain%svertex(nv)%Veloc(1)
              write (63,*) nv+1, Tdomain%svertex(nv)%Veloc(2)
              count(nv) = 1
           endif
        enddo
        !  endif
     enddo


     close(61)
     close(62)
     close(63)

     deallocate(count)
  endif !Order Veloc

  if (Tdomain%Field_Order(2)) then !Stress Order

     call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
     if (.not. (allocated(NbElem))) allocate(NbElem(0:nbprocs-1)) !stocking the number of non PML elements
     call MPI_ALLGATHER(Tdomain%n_elem_nonPML,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
          MPI_COMM_WORLD,code)

     if (.not. allocated(NbGlob)) allocate(NbGlob(0:nbprocs-1))
     call MPI_ALLGATHER(Tdomain%n_glob_points,1,MPI_INTEGER,NbGlob,1,MPI_INTEGER, &
          MPI_COMM_WORLD,code)
!!$     if (it==Tdomain%TimeD%NtimeMax-1) print *,'STRESS',NbElem(rg),NbGlob(rg),rg
     if (icount==1) print *,'STRESS',NbElem(rg),NbGlob(rg),rg
     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,code)
     call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,code)

     write (fnamef,"(a,I4.4,a)") "SField/StressField",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,&
          MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

     SZG=3*Tdomain%n_glob_points

     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_elem_nonPML,1,MPI_INTEGER,statut,code) !mpi_write2
     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_glob_points,1,MPI_INTEGER,statut,code) !mpi_write3
     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%GlobCoord,SZG,MPI_DOUBLE_PRECISION,statut,code) !mpi_write4
!!$if (it==Tdomain%TimeD%NtimeMax-1) print *, 'SIGMA',size(Tdomain%GlobCoord),SZG,rg
     PosFile = (2 * nbprocs) * NbOctInt !Jump over mpi_write1-2-3
     do i1 = 0, nbprocs-1 !Jump over  mpi_write4
        PosFile = PosFile + 3*NbGlob(i1)*NbOctReal
     enddo

     do n = 0, Tdomain%n_elem-1

        if (.not. Tdomain%specel(n)%PML) then
           ngllx = Tdomain%specel(n)%ngllx
           nglly = Tdomain%specel(n)%nglly
           ngllz = Tdomain%specel(n)%ngllz
           mat = Tdomain%specel(n)%mat_index
           if (.not. allocated(GlobId)) allocate(GlobId(0:ngllx-1,0:nglly-1,0:ngllz-1))
           GlobId=Tdomain%specel(n)%Iglobnum

           call Displacement_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
           call Displacement_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
           call Displacement_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)

           if (.not. allocated (dUx_dxi)) allocate (dUx_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_deta)) allocate (dUx_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUx_dzeta)) allocate (dUx_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dxi)) allocate (dUy_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_deta)) allocate (dUy_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUy_dzeta)) allocate (dUy_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dxi)) allocate (dUz_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_deta)) allocate (dUz_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
           if (.not. allocated (dUz_dzeta)) allocate (dUz_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))

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

           allocate (xix(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (xiy(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (xiz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etax(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etay(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (etaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (zetax(0:ngllx-1,0:nglly-1,0:ngllz-1))
           allocate (zetay(0:ngllx-1,0:nglly-1,0:ngllz-1)) 
           allocate (zetaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
!if (rg==3) print *,'rg3--', Tdomain%specel(30)%InvGrad (0,0,:, 0,0)
           xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)

           xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
           xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)

           etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
           etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
           etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

           zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
           zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
           zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)


           if (.not. allocated (stress_sigma)) allocate (stress_sigma(0:ngllx-1,0:nglly-1,0:ngllz-1,1:6))           
           if (.not. allocated (strain)) allocate (strain(0:ngllx-1,0:nglly-1,0:ngllz-1,1:6))           
           strain(:,:,:,1)=&
                dUx_dxi*xix+dUx_deta*etax+dUx_dzeta*zetax
           strain(:,:,:,2)=&
                dUy_dxi*xix+dUy_deta*etay+dUy_dzeta*zetay
           strain(:,:,:,3)=&
                dUz_dxi*xiz+dUz_deta*etaz+dUz_dzeta*zetaz
           strain(:,:,:,4)=&
                (dUx_dxi*xiy+dUx_deta*etay+dUx_dzeta*zetay&
                +dUy_dxi*xix+dUy_deta*etax+dUy_dzeta*zetax)/2.
           strain(:,:,:,5)=&
                (dUx_dxi*xiz+dUx_deta*etaz+dUx_dzeta*zetaz&
                +dUz_dxi*xix+dUz_deta*etax+dUz_dzeta*zetax)/2.
           strain(:,:,:,6)=&
                (dUy_dxi*xiz+dUy_deta*etaz+dUy_dzeta*zetaz&
                +dUz_dxi*xiy+dUz_deta*etay+dUz_dzeta*zetay)/2.

           deallocate (xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz)

           do n_x=0,ngllx-1
              do n_y=0,nglly-1
                 do n_z=0,ngllz-1
                    call DGEMV('N', 6, 6, 1.d0,Tdomain%specel(n)%c(:,:,n_x,n_y,n_z), &
                    6,strain(n_x,n_y,n_z,:) , 1, 0.d0, stress_sigma(n_x,n_y,n_z,:), 1)
                   !DGEMV(TRANSA, M, N, ALPHA,                   A,
                   !!LDA,   X,             INCX, BETA,         Y,                  INCY)
                 enddo
              enddo
           enddo

           SZG=ngllx*nglly*ngllz
           SZG6=6*SZG
           PosFileR = PosFile
           PosFileI = PosFile
           do i1 = 0, rg-1 !Previous rg
              PosFileI = PosFileI+NbElem(i1)*NbOctInt*SZG
              PosFileR = PosFileR+NbElem(i1)*NbOctReal*SZG6
           enddo


           do i1 = 0, nbprocs-1
              PosFileR = PosFileR + NbElem(i1)*NbOctInt*SZG ! C2write "waits for" GlobId
           end do
           PosFileR = PosFileR + n*NbOctReal*SZG*6 !stress_sigma of previous elts
           PosFileI = PosFileI + n*NbOctInt*SZG !GlobId of previous elts
           call MPI_FILE_WRITE_AT(desc,PosFileI,GlobId,SZG,MPI_INTEGER,statut,code)
           call MPI_FILE_WRITE_AT(desc,PosFileR,stress_sigma(:,:,:,:),SZG6,MPI_DOUBLE_PRECISION,statut,code)
        endif
     enddo
     if ((rg==nbprocs-1)) then
        PosFile=PosFileR+SZG6*NbOctReal
        call MPI_FILE_WRITE_AT(desc,PosFile,nbprocs,1,MPI_INTEGER,statut,code)  !mpi_write1
     endif
     call MPI_FILE_CLOSE(desc,code)

  endif !Order Stress



  if (Tdomain%Field_Order(3)) then !Energy Order




  endif !Order Energy

  if (Tdomain%Field_Order(4)) then !Traveltime Order
     do n = 0, Tdomain%n_elem-1
        if (.not. Tdomain%specel(n)%PML) then
           ngllx = Tdomain%specel(n)%ngllx
           nglly = Tdomain%specel(n)%nglly
           ngllz = Tdomain%specel(n)%ngllz
           if (.not. allocated(GlobId)) allocate(GlobId(0:ngllx-1,0:nglly-1,0:ngllz-1))
           GlobId=Tdomain%specel(n)%Iglobnum
           call Velocity_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
           call Velocity_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
           call Velocity_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)
           if (.not. Tdomain%replay) then
              do n_z = 0, ngllz -1 
                 do n_y = 0, nglly -1 
                    do n_x = 0, ngllx -1
                       AbsVeloc=sqrt(Tdomain%specel(n)%Veloc(n_x,n_y,n_z,0)**2+&
                            Tdomain%specel(n)%Veloc(n_x,n_y,n_z,1)**2+&
                            Tdomain%specel(n)%Veloc(n_x,n_y,n_z,2)**2)
                       if (AbsVeloc .gt. Tdomain%specel(n)%MaxAbsVeloc (n_x,n_y,n_z)) &
                       Tdomain%specel(n)%MaxAbsVeloc (n_x,n_y,n_z)=AbsVeloc
                    enddo
                 enddo
              enddo
           else
              do n_z = 0, ngllz -1 
                 do n_y = 0, nglly -1 
                    do n_x = 0, ngllx -1
                       AbsVeloc=sqrt(Tdomain%specel(n)%Veloc(n_x,n_y,n_z,0)**2+&
                            Tdomain%specel(n)%Veloc(n_x,n_y,n_z,1)**2+&
                            Tdomain%specel(n)%Veloc(n_x,n_y,n_z,2)**2)
                       do l=1,size(Tdomain%Ratios)
                          if (.not. Tdomain%specel(n)%TravelTimeFound(n_x,n_y,n_z,l)) then
                             if (.not. (AbsVeloc .lt. (Tdomain%Ratios(l)*Tdomain%specel(n)%MaxAbsVeloc (n_x, n_y, n_z)))) then
!!!                                Tdomain%specel(n)%TravelTime(n_x,n_y,n_z,l)=it*Tdomain%sSubdomain(Tdomain%specel(n)%mat_index)%Dt
                                Tdomain%specel(n)%TravelTimeFound(n_x,n_y,n_z,l)=.true.
!!!                                print *, 'TravelTime(nx ny nz rat)max it dt ',Tdomain%specel(n)%TravelTime(n_x,n_y,n_z,l),&
!!!                                Tdomain%specel(n)%MaxAbsVeloc (n_x, n_y, n_z),it,Tdomain%sSubdomain(Tdomain%specel(n)%mat_index)%Dt
                             endif
                          endif
                       enddo
!!!                       if (it==0) Tdomain%specel(n)%TravelTime(n_x,n_y,n_z,size(Tdomain%Ratios)+1)=Tdomain%specel(n)%MaxAbsVeloc (n_x,n_y,n_z)
                    enddo
                 enddo
              enddo
           endif !if replay

        endif
     enddo

     !Writting
     if ((it==(Tdomain%TimeD%NtimeMax-1)) .and. (Tdomain%replay)) then

        call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
        allocate(NbElem(0:nbprocs-1)) !stocking the number of non PML elements
        call MPI_ALLGATHER(Tdomain%n_elem_nonPML,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
             MPI_COMM_WORLD,code)

        allocate(NbGlob(0:nbprocs-1))
        call MPI_ALLGATHER(Tdomain%n_glob_points,1,MPI_INTEGER,NbGlob,1,MPI_INTEGER, &
             MPI_COMM_WORLD,code)
        print *,'TIME1',NbElem(rg),NbGlob(rg),rg
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,code)
        call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,code)

        write (fnamef,"(a)") "SField/TravelTimeField.out"     
        call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,&
             MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

        SZG=3*Tdomain%n_glob_points

        call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_elem_nonPML,1,MPI_INTEGER,statut,code) !mpi_write2
        call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_glob_points,1,MPI_INTEGER,statut,code) !mpi_write3
        call MPI_FILE_WRITE_ORDERED(desc,Tdomain%GlobCoord,SZG,MPI_DOUBLE_PRECISION,statut,code) !mpi_write4
        print *,'TIME2',size(Tdomain%GlobCoord),SZG,rg
        PosFile = (2 * nbprocs) * NbOctInt !Jump over mpi_write1-2-3
        do i1 = 0, nbprocs-1 !Jump over  mpi_write4
           PosFile = PosFile + 3*NbGlob(i1)*NbOctReal
        enddo

        do n = 0, Tdomain%n_elem-1
           if (.not. Tdomain%specel(n)%PML) then
              ngllx = Tdomain%specel(n)%ngllx
              nglly = Tdomain%specel(n)%nglly
              ngllz = Tdomain%specel(n)%ngllz
              GlobId=Tdomain%specel(n)%Iglobnum
              SZG=ngllx*nglly*ngllz
              SZG6=(size(Tdomain%Ratios)+1)*SZG
              PosFileR = PosFile
              PosFileI = PosFile
              do i1 = 0, rg-1 !Previous rg
                 PosFileI = PosFileI+NbElem(i1)*NbOctInt*SZG
                 PosFileR = PosFileR+NbElem(i1)*NbOctReal*SZG6
              enddo
              do i1 = 0, nbprocs-1
                 PosFileR = PosFileR + NbElem(i1)*NbOctInt*SZG ! C2write "waits for" GlobId
              end do
              PosFileR = PosFileR + n*NbOctReal*SZG6 ! Traveltime of previous elts
              PosFileI = PosFileI + n*NbOctInt*SZG !GlobId of previous elts
              call MPI_FILE_WRITE_AT(desc,PosFileI,GlobId,SZG,MPI_INTEGER,statut,code)
!!!              call MPI_FILE_WRITE_AT(desc,PosFileR,Tdomain%specel(n)%Traveltime(:,:,:,:),&
!!!              SZG6,MPI_DOUBLE_PRECISION,statut,code)
           endif
        enddo
        if ((rg==nbprocs-1)) then
           PosFile=PosFileR+SZG6*NbOctReal
           print *, 'TIME',Posfile       
           call MPI_FILE_WRITE_AT(desc,PosFile,nbprocs,1,MPI_INTEGER,statut,code)  !mpi_write1
        endif
        call MPI_FILE_CLOSE(desc,code)

     endif


  endif !Order Traveltime

  if (Tdomain%Field_Order(5)) then !Veloc Order at all of GLL points
     call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
     if (.not. (allocated(NbElem))) allocate(NbElem(0:nbprocs-1)) !stocking the number of non PML elements
     call MPI_ALLGATHER(Tdomain%n_elem_nonPML,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
          MPI_COMM_WORLD,code)
!!$call MPI_ALLGATHER(Tdomain%n_elem,1,MPI_INTEGER,NbElem,1,MPI_INTEGER, &
!!$          MPI_COMM_WORLD,code)
     if (.not. allocated(NbGlob)) allocate(NbGlob(0:nbprocs-1))
     call MPI_ALLGATHER(Tdomain%n_glob_points,1,MPI_INTEGER,NbGlob,1,MPI_INTEGER, &
          MPI_COMM_WORLD,code)
!!$     if (it==Tdomain%TimeD%NtimeMax-1) print *,'STRESS',NbElem(rg),NbGlob(rg),rg
     if (icount==1) print *,'STRESS',NbElem(rg),NbGlob(rg),rg
     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,code)
     call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,code)

     write (fnamef,"(a,I4.4,a)") "SField/VelocGLLsField",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,&
          MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

     SZG=3*Tdomain%n_glob_points

     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_elem_nonPML,1,MPI_INTEGER,statut,code) !mpi_write2
!!$ call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_elem,1,MPI_INTEGER,statut,code) !mpi_write2
     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%n_glob_points,1,MPI_INTEGER,statut,code) !mpi_write3
     call MPI_FILE_WRITE_ORDERED(desc,Tdomain%GlobCoord,SZG,MPI_DOUBLE_PRECISION,statut,code) !mpi_write4
!!$if (it==Tdomain%TimeD%NtimeMax-1) print *, 'SIGMA',size(Tdomain%GlobCoord),SZG,rg
     PosFile = (2 * nbprocs) * NbOctInt !Jump over mpi_write1-2-3
     do i1 = 0, nbprocs-1 !Jump over  mpi_write4
        PosFile = PosFile + 3*NbGlob(i1)*NbOctReal
     enddo
     do n = 0, Tdomain%n_elem-1

        if (.not. Tdomain%specel(n)%PML) then
           ngllx = Tdomain%specel(n)%ngllx
           nglly = Tdomain%specel(n)%nglly
           ngllz = Tdomain%specel(n)%ngllz
           mat = Tdomain%specel(n)%mat_index
           if (.not. allocated(GlobId)) allocate(GlobId(0:ngllx-1,0:nglly-1,0:ngllz-1))
           GlobId=Tdomain%specel(n)%Iglobnum

           call Velocity_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
           call Velocity_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
           call Velocity_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)

           if (.not. allocated (Veloc_GLLs)) allocate (Veloc_GLLs(0:ngllx-1,0:nglly-1,0:ngllz-1,1:4))  
!!$       ! Veloc_GLLs(:,:,:,1:3)=Tdomain%specel(n)%Veloc
           do n_x=0,ngllx-1
              do n_y=0,nglly-1
                 do n_z=0,ngllz-1
                    Veloc_GLLs(n_x,n_y,n_z,1)=Tdomain%specel(n)%Veloc(n_x,n_y,n_z,0)
                    Veloc_GLLs(n_x,n_y,n_z,2)=Tdomain%specel(n)%Veloc(n_x,n_y,n_z,1)
                    Veloc_GLLs(n_x,n_y,n_z,3)=Tdomain%specel(n)%Veloc(n_x,n_y,n_z,2)
                    Veloc_GLLs(n_x,n_y,n_z,4)=sqrt(Veloc_GLLs(n_x,n_y,n_z,1)**2&
                    +Veloc_GLLs(n_x,n_y,n_z,2)**2+Veloc_GLLs(n_x,n_y,n_z,3)**2)
!!$                   ! Veloc_GLLs(n_x,n_y,n_z,1:6)=0.
                 enddo
              enddo
           enddo
           SZG=ngllx*nglly*ngllz
           SZG6=4*SZG
           PosFileR = PosFile
           PosFileI = PosFile
           do i1 = 0, rg-1 !Previous rg
              PosFileI = PosFileI+NbElem(i1)*NbOctInt*SZG
              PosFileR = PosFileR+NbElem(i1)*NbOctReal*SZG6
           enddo


           do i1 = 0, nbprocs-1
              PosFileR = PosFileR + NbElem(i1)*NbOctInt*SZG ! C2write "waits for" GlobId
           end do
           PosFileR = PosFileR + n*NbOctReal*SZG*4 !stress_sigma of previous elts
           PosFileI = PosFileI + n*NbOctInt*SZG !GlobId of previous elts
           call MPI_FILE_WRITE_AT(desc,PosFileI,GlobId,SZG,MPI_INTEGER,statut,code)
           call MPI_FILE_WRITE_AT(desc,PosFileR,Veloc_GLLs(:,:,:,:),SZG6,MPI_DOUBLE_PRECISION,statut,code)

        endif
     enddo
     if ((rg==nbprocs-1)) then
        PosFile=PosFileR+SZG6*NbOctReal
        if (it==Tdomain%TimeD%NtimeMax-1) print *,'VelocGLLs',PosFile
        call MPI_FILE_WRITE_AT(desc,PosFile,nbprocs,1,MPI_INTEGER,statut,code)  !mpi_write1
     endif
     call MPI_FILE_CLOSE(desc,code)
  endif

  if (Tdomain%Field_Order(6)) then
     call save_slices(Tdomain,it,rg,icount)
  endif
  if (Tdomain%Field_Order(7)) then
     call save_slices_enerPS(Tdomain,it,rg,icount)
  endif



  return

end subroutine savefield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from vertex to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displacement_Collection_From_Vertex(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)


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
end subroutine Displacement_Collection_From_Vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Face to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displacement_Collection_From_Face(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)
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
end subroutine Displacement_Collection_From_Face
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Edges to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Displacement_Collection_From_Edge(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)

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
end subroutine Displacement_Collection_From_Edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from vertex to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Velocity_Collection_From_Vertex(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)


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
end subroutine Velocity_Collection_From_Vertex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Face to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Velocity_Collection_From_Face(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)
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
end subroutine Velocity_Collection_From_Face
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SUB-SUBROUTINE collecting the displacement from Edges to element !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine Velocity_Collection_From_Edge(s_Tdomain,s_n,s_ngllx,s_nglly,s_ngllz)

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
end subroutine Velocity_Collection_From_Edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
