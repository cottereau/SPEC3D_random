subroutine save_slices_enerPS(Tdomain,it,rg,icount)
  !Veloc Order at all of GLL points of 3 slices in x-,y-,z- directions 


  use sdomain

  implicit none

  include 'mpif.h'


  type (domain), intent (INOUT):: Tdomain
  integer, intent (IN) :: it,rg,icount

  ! local variables
  integer :: i,n,nv,nbvert,kcount,nv_aus
  integer, dimension(:), allocatable:: count
  character (len=100) :: fnamef

  integer, dimension(:),allocatable ::NbElem_slice0,NbElem_slice1,NbElem_slice2, NbGlob 
  integer::SZG,SZG3,n_x,n_y,n_z,i1,NbOctInt, NbOctReal,desc,code,nbprocs,ngllx,nglly,ngllz,mat,l,nx1,nx2,ipoint,ns,ne,descGeo,ne_above,ne_left


  integer,dimension(:,:,:),allocatable::GlobId
!!$  real, dimension (:,:,:,:), allocatable :: stress_sigma,strain,Veloc_GLLs
  real, dimension (:,:,:), allocatable :: EnerGLL_elem_Slice,GeoGLL_elem_Slice
  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer (kind=MPI_OFFSET_KIND) :: PosFile, PosFileI, PosFileVelR,PosFileGeoR,seek_offset=0
  integer,dimension(:,:),pointer::NbElem_inSlices

  real, dimension(:,:,:), allocatable :: dVx_dxi, dVx_deta, dVx_dzeta, &
       dVy_dxi, dVy_deta, dVy_dzeta, &
       dVz_dxi, dVz_deta, dVz_dzeta, &
       dVx_dx, dVx_dy, dVx_dz, &
       dVy_dx, dVy_dy, dVy_dz, &
       dVz_dx, dVz_dy, dVz_dz, &
       xix ,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz
  real, dimension(:,:,:,:), allocatable :: Ediv_curl 



  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbprocs,code)
  call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,NbOctReal,code)
  call MPI_TYPE_SIZE(MPI_INTEGER,NbOctInt,code)


  !print *,NbOctInt,NbOctReal

  PosFileVelR=0
  PosFileGeoR=0

  allocate(NbElem_inSlices(Tdomain%nbSlices,0:nbprocs-1))

  NbElem_inSlices=0


  if (icount==1) then


     write (fnamef,"(a)") "SField/GeoSlices.out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,descGeo,code)
     write (fnamef,"(a,I4.4,a)") "SField/EnerPSSlices",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)
	print *, 'PScheck0', Tdomain%nbSlices
     do ns=1,Tdomain%nbSlices
       print *,'PScheck1', Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize),rg 
        call MPI_ALLGATHER(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize),1,MPI_INTEGER,NbElem_inSlices(ns,:),1,MPI_INTEGER, MPI_COMM_WORLD,code)
        if (ns .gt. 1 ) then 
           ne_above=sum(NbElem_inSlices(1:ns-1,:))
        else
           ne_above=0
        endif
        
        if (.not.(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize) .lt. 1)) then
               
           do ne=1,Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize)
              if ((rg .gt. 0)) then 
                 ne_left=sum(NbElem_inSlices(ns,0:rg-1))+(ne-1)
              else
                 ne_left=ne-1
              endif
              n=Tdomain%I_elem_inSlices(ns,ne)

              ngllx=Tdomain%specel(n)%ngllx
              nglly=Tdomain%specel(n)%nglly
              ngllz=Tdomain%specel(n)%ngllz


	!print *, ngllx, nglly,ngllz
              call Velocity_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
              call Velocity_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
              call Velocity_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)

mat = Tdomain%specel(n)%mat_index

              if (.not. allocated (dVx_dxi)) allocate (dVx_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_deta)) allocate (dVx_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dzeta)) allocate (dVx_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dxi)) allocate (dVy_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_deta)) allocate (dVy_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dzeta)) allocate (dVy_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dxi)) allocate (dVz_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_deta)) allocate (dVz_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dzeta)) allocate (dVz_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))


              if (.not. allocated (dVx_dx)) allocate (dVx_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dy)) allocate (dVx_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dz)) allocate (dVx_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dx)) allocate (dVy_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dy)) allocate (dVy_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dz)) allocate (dVy_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dx)) allocate (dVz_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dy)) allocate (dVz_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dz)) allocate (dVz_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))

  


              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                   ngllx, Tdomain%specel(n)%Veloc(0,0,0,0), ngllx, 0., dVx_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,0),   &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVx_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,0),  &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVx_dzeta, ngllx*nglly)
              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                   ngllx,Tdomain%specel(n)%Veloc(0,0,0,1), ngllx, 0., dVy_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,1),   &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVy_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,1),   &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVy_dzeta, ngllx*nglly)
              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,    &
                   ngllx,Tdomain%specel(n)%Veloc(0,0,0,2), ngllx, 0., dVz_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,2),    &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVz_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,2),   &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVz_dzeta, ngllx*nglly)

            
              allocate (xix(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (xiy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (xiz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etax(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etay(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (zetax(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (zetay(0:ngllx-1,0:nglly-1,0:ngllz-1)) 
              allocate (zetaz(0:ngllx-1,0:nglly-1,0:ngllz-1))


           xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)

           xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
           xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)

           etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
           etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
           etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)


           zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
           zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
           zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)


              dVx_dx=dVx_dxi*xix+dVx_deta*etax+dVx_dzeta*zetax
              dVy_dx=dVy_dxi*xix+dVy_deta*etax+dVy_dzeta*zetax
              dVz_dx=dVz_dxi*xix+dVz_deta*etax+dVz_dzeta*zetax

              dVx_dy=dVx_dxi*xiy+dVx_deta*etay+dVx_dzeta*zetay
              dVy_dy=dVy_dxi*xiy+dVy_deta*etay+dVy_dzeta*zetay
              dVz_dy=dVz_dxi*xiy+dVz_deta*etay+dVz_dzeta*zetay

              dVx_dz=dVx_dxi*xiz+dVx_deta*etaz+dVx_dzeta*zetaz
              dVy_dz=dVy_dxi*xiz+dVy_deta*etaz+dVy_dzeta*zetaz
              dVz_dz=dVz_dxi*xiz+dVz_deta*etaz+dVz_dzeta*zetaz
              deallocate(dVx_dxi,dVx_deta,dVx_dzeta,&
                   dVy_dxi,dVy_deta,dVy_dzeta,&
                   dVz_dxi,dVz_deta,dVz_dzeta,&
                   xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz)
              
allocate(Ediv_curl(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
Ediv_curl(:,:,:,0)=(dVx_dx+dVy_dy+dVz_dz)**2
Ediv_curl(:,:,:,1)=(dVy_dz-dVz_dy)**2+(dVz_dx-dVx_dz)**2+(dVx_dy-dVy_dx)**2
!!$              allocate (Ediv(0:ngllx-1,0:nglly-1,0:ngllz-1))
!!$              allocate (Ecurl(0:ngllx-1,0:nglly-1,0:ngllz-1))
!!$              Ediv(:,:,:)=(dVx_dx+dVy_dy+dVz_dz)**2
!!$              Ediv(:,:,:)=(dVy_dz-dVz_dy)**2+(dVz_dx-dVx_dz)**2+(dVx_dy-dVy_dx)**2
              deallocate(dVx_dx,dVx_dy,dVx_dz,&
                   dVy_dx,dVy_dy,dVy_dz,&
                   dVz_dx,dVz_dy,dVz_dz)
              
              select case(int(Tdomain%Slices(ns,1)))
              case (0)
                 allocate(GeoGLL_elem_Slice(0:nglly-1,0:ngllz-1,2))
                 allocate(EnerGLL_elem_Slice(0:nglly-1,0:ngllz-1,0:1))
                 SZG=nglly*ngllz
                 do nx1=0,nglly-1
                    do nx2=0,ngllz-1
                       ipoint = Tdomain%specel(n)%Iglobnum(0,nx1,nx2)

                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (1,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (2,ipoint) 
                       EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(0,nx1,nx2,:)
                      
!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(0,nx1,nx2)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(0,nx1,nx2)
                    enddo
                 enddo

              case (1)
                 allocate(GeoGLL_elem_Slice(0:ngllx-1,0:ngllz-1,2))
                 allocate(EnerGLL_elem_Slice(0:ngllx-1,0:ngllz-1,0:1))
                 SZG=ngllx*ngllz
                 do nx1=0,ngllx-1
                    do nx2=0,ngllz-1
                       ipoint = Tdomain%specel(n)%Iglobnum(nx1,0,nx2)

                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (0,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (2,ipoint)
                       EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(nx1,0,nx2,:)
!!$
!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(nx1,0,nx2)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(nx1,0,nx2)

                    enddo
                 enddo

              case (2)

                 allocate(GeoGLL_elem_Slice(0:ngllx-1,0:nglly-1,2))
                 allocate(EnerGLL_elem_Slice(0:ngllx-1,0:nglly-1,0:1))
                 SZG=nglly*ngllx
                 do nx1=0,ngllx-1
                    do nx2=0,nglly-1
                       ipoint = Tdomain%specel(n)%Iglobnum(nx1,nx2,ngllz-1)
                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (0,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (1,ipoint) 
                       EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(nx1,nx2,ngllz-1,:)

!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(nx1,nx2,ngllz-1)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(nx1,nx2,ngllz-1)
                    enddo
                 enddo
              endselect
              PosFileGeoR=(ne_above+ne_left)*NbOctReal*SZG*2
              PosFileVelR=(ne_above+ne_left)*NbOctReal*SZG*2
              call MPI_FILE_WRITE_AT(descGeo,PosFileGeoR,GeoGLL_elem_Slice(:,:,:),2*SZG,MPI_DOUBLE_PRECISION,statut,code)
              call MPI_FILE_WRITE_AT(desc,PosFileVelR,EnerGLL_elem_Slice(:,:,:),2*SZG,MPI_DOUBLE_PRECISION,statut,code)           
              if (Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize).gt.0) deallocate(GeoGLL_elem_Slice,EnerGLL_elem_Slice,Ediv_curl)
           enddo
        endif
     enddo
     call MPI_BARRIER(MPI_COMM_WORLD,code)
     if ((rg==0)) then
        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_END,code)
        !call system('ls -lstr SField/GeoSlices.out')
        SZG=Tdomain%nbSlices*nbprocs

        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo,Nbelem_inSlices(:,:),SZG,MPI_INTEGER,statut,code)
        !call system('ls -lstr SField/GeoSlices.out')
        !----- tic-toc-step
        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo, Tdomain%TimeD%time_snapshots_tic,1,MPI_DOUBLE_PRECISION,statut,code)
        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo, Tdomain%TimeD%time_snapshots_step,1,MPI_DOUBLE_PRECISION,statut,code)
        !----

        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo,int(Tdomain%Slices(:,1)),Tdomain%nbSlices,MPI_INTEGER,statut,code)
        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo,Tdomain%Slices(:,2),Tdomain%nbSlices,MPI_DOUBLE_PRECISION,statut,code)
        !call system('ls -lstr SField/GeoSlices.out') 
        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo,nbprocs,1,MPI_INTEGER,statut,code) 

        call MPI_FILE_SEEK(descGeo,seek_offset,MPI_SEEK_CUR,code)
        call MPI_FILE_WRITE(descGeo,Tdomain%nbSlices,1,MPI_INTEGER,statut,code) 

     endif
     call MPI_FILE_CLOSE(descGeo,code)
     call MPI_FILE_CLOSE(desc,code)


  else
     write (fnamef,"(a,I4.4,a)") "SField/EnerPSSlices",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

     do ns=1,Tdomain%nbSlices

        call MPI_ALLGATHER(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize),1,MPI_INTEGER,NbElem_inSlices(ns,:),1,MPI_INTEGER, MPI_COMM_WORLD,code)
        if (ns .gt. 1 ) then 
           ne_above=sum(NbElem_inSlices(1:ns-1,:))
        else
           ne_above=0
        endif

        if (.not.(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize) .lt. 1)) then
           do ne=1,Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize)
              if ((rg .gt. 0)) then 
                 ne_left=sum(NbElem_inSlices(ns,0:rg-1))+(ne-1)
              else
                 ne_left=ne-1
              endif
              n=Tdomain%I_elem_inSlices(ns,ne)
              ngllx=Tdomain%specel(n)%ngllx
              nglly=Tdomain%specel(n)%nglly
              ngllz=Tdomain%specel(n)%ngllz
              !print *, ngllx, nglly,ngllz
              call Velocity_Collection_From_Vertex(Tdomain,n,ngllx,nglly,ngllz)
              call Velocity_Collection_From_Edge(Tdomain,n,ngllx,nglly,ngllz)
              call Velocity_Collection_From_Face(Tdomain,n,ngllx,nglly,ngllz)
  mat = Tdomain%specel(n)%mat_index

              if (.not. allocated (dVx_dxi)) allocate (dVx_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_deta)) allocate (dVx_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dzeta)) allocate (dVx_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dxi)) allocate (dVy_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_deta)) allocate (dVy_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dzeta)) allocate (dVy_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dxi)) allocate (dVz_dxi(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_deta)) allocate (dVz_deta(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dzeta)) allocate (dVz_dzeta(0:ngllx-1,0:nglly-1,0:ngllz-1))


              if (.not. allocated (dVx_dx)) allocate (dVx_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dy)) allocate (dVx_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVx_dz)) allocate (dVx_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dx)) allocate (dVy_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dy)) allocate (dVy_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVy_dz)) allocate (dVy_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dx)) allocate (dVz_dx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dy)) allocate (dVz_dy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              if (.not. allocated (dVz_dz)) allocate (dVz_dz(0:ngllx-1,0:nglly-1,0:ngllz-1))



              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                   ngllx, Tdomain%specel(n)%Veloc(0,0,0,0), ngllx, 0., dVx_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,0),   &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVx_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,0),  &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVx_dzeta, ngllx*nglly)
              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,   &
                   ngllx,Tdomain%specel(n)%Veloc(0,0,0,1), ngllx, 0., dVy_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,1),   &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVy_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,1),   &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVy_dzeta, ngllx*nglly)
              call DGEMM ( 'N', 'N', ngllx, nglly*ngllz, ngllx, 1., Tdomain%sSubdomain(mat)%hTprimex,    &
                   ngllx,Tdomain%specel(n)%Veloc(0,0,0,2), ngllx, 0., dVz_dxi, ngllx)
              do n_z = 0,ngllz-1
                 call DGEMM ( 'N', 'N', ngllx, nglly, nglly, 1., Tdomain%specel(n)%Veloc(0,0,n_z,2),    &
                      ngllx, Tdomain%sSubdomain(mat)%hprimey, nglly, 0., dVz_deta(0,0,n_z), ngllx)
              enddo
              call DGEMM ( 'N', 'N', ngllx*nglly, ngllz, ngllz, 1., Tdomain%specel(n)%Veloc(0,0,0,2),   &
                   ngllx*nglly, Tdomain%sSubdomain(mat)%hprimez, ngllz, 0., dVz_dzeta, ngllx*nglly)

              allocate (xix(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (xiy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (xiz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etax(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etay(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (etaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (zetax(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate (zetay(0:ngllx-1,0:nglly-1,0:ngllz-1)) 
              allocate (zetaz(0:ngllx-1,0:nglly-1,0:ngllz-1))


              xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)

              xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
              xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)

              etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
              etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
              etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

              zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
              zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
              zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)


              dVx_dx=dVx_dxi*xix+dVx_deta*etax+dVx_dzeta*zetax
              dVy_dx=dVy_dxi*xix+dVy_deta*etax+dVy_dzeta*zetax
              dVz_dx=dVz_dxi*xix+dVz_deta*etax+dVz_dzeta*zetax

              dVx_dy=dVx_dxi*xiy+dVx_deta*etay+dVx_dzeta*zetay
              dVy_dy=dVy_dxi*xiy+dVy_deta*etay+dVy_dzeta*zetay
              dVz_dy=dVz_dxi*xiy+dVz_deta*etay+dVz_dzeta*zetay

              dVx_dz=dVx_dxi*xiz+dVx_deta*etaz+dVx_dzeta*zetaz
              dVy_dz=dVy_dxi*xiz+dVy_deta*etaz+dVy_dzeta*zetaz
              dVz_dz=dVz_dxi*xiz+dVz_deta*etaz+dVz_dzeta*zetaz

              allocate (Ediv_curl(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
              Ediv_curl(:,:,:,0)=(dVx_dx+dVy_dy+dVz_dz)**2
              Ediv_curl(:,:,:,1)=(dVy_dz-dVz_dy)**2+(dVz_dx-dVx_dz)**2+(dVx_dy-dVy_dx)**2
              deallocate(dVx_dxi,dVx_deta,dVx_dzeta,&
                   dVy_dxi,dVy_deta,dVy_dzeta,&
                   dVz_dxi,dVz_deta,dVz_dzeta,&
                   xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz)

!!$              allocate (Ediv(0:ngllx-1,0:nglly-1,0:ngllz-1))
!!$              allocate (Ecurl(0:ngllx-1,0:nglly-1,0:ngllz-1))
!!$
!!$              Ediv(:,:,:)=(dVx_dx+dVy_dy+dVz_dz)**2
!!$              Ecurl(:,:,:)=(dVy_dz-dVz_dy)**2+(dVz_dx-dVx_dz)**2+(dVx_dy-dVy_dx)**2
!!$



              select case(int(Tdomain%Slices(ns,1)))
              case (0)

                 allocate(EnerGLL_elem_Slice(0:nglly-1,0:ngllz-1,0:2))
                 SZG=nglly*ngllz
                 do nx1=0,nglly-1
                    do nx2=0,ngllz-1

                    EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(0,nx1,nx2,:)
!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(0,nx1,nx2)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(0,nx1,nx2)
                       !print *,'slices_check1',SZG
                    enddo
                 enddo

              case (1)

                 allocate(EnerGLL_elem_Slice(0:ngllx-1,0:ngllz-1,0:2))
                 SZG=ngllx*ngllz
                 do nx1=0,ngllx-1
                    do nx2=0,ngllz-1

                    EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(nx1,0,nx2,:)
!!$
!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(nx1,0,nx2)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(nx1,0,nx2)
                    enddo
                 enddo

              case (2)

                 allocate(EnerGLL_elem_Slice(0:ngllx-1,0:nglly-1,0:2))
                 SZG=nglly*ngllx
                 do nx1=0,ngllx-1
                    do nx2=0,nglly-1

                    EnerGLL_elem_Slice(nx1,nx2,:)=Ediv_curl(nx1,nx2,ngllz-1,:)
!!$
!!$                       EnerGLL_elem_Slice(nx1,nx2,0)=Ediv(nx1,nx2,ngllz-1)
!!$                       EnerGLL_elem_Slice(nx1,nx2,1)=Ecurl(nx1,nx2,ngllz-1)
                    enddo
                 enddo
              endselect

              PosFileVelR=(ne_above+ne_left)*NbOctReal*SZG*2
              call MPI_FILE_WRITE_AT(desc,PosFileVelR,EnerGLL_elem_Slice(:,:,:),2*SZG,MPI_DOUBLE_PRECISION,statut,code)           
              if (Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize).gt.0) deallocate(EnerGLL_elem_Slice,Ediv_curl)
           enddo
        endif
     enddo
     call MPI_BARRIER(MPI_COMM_WORLD,code)
     call MPI_FILE_CLOSE(desc,code)

  endif

  deallocate(NbElem_inSlices)

  return

end subroutine save_slices_enerPS


