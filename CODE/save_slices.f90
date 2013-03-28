subroutine save_slices(Tdomain,it,rg,icount)
  !Veloc Order at all of GLL points of 3 slices in x-,y-,z- directions 


  use sdomain

  implicit none

  include 'mpif.h'


  type (domain), intent (INOUT):: Tdomain
  integer, intent (IN) :: it,rg,icount

  ! local variables
  integer :: i,n,nv,nbvert,kcount,nv_aus,ne_above,ne_left
  integer, dimension(:), allocatable:: count
  character (len=100) :: fnamef

  integer, dimension(:),allocatable ::NbElem_slice0,NbElem_slice1,NbElem_slice2, NbGlob 
  integer::SZG,SZG3,n_x,n_y,n_z,i1,NbOctInt, NbOctReal,desc,code,nbprocs,ngllx,nglly,ngllz,mat,l,nx1,nx2,ipoint,ns,ne,descGeo,ne_above,ne_left


  integer,dimension(:,:,:),allocatable::GlobId
!!$  real, dimension (:,:,:,:), allocatable :: stress_sigma,strain,Veloc_GLLs
  real, dimension (:,:,:), allocatable :: VelocGLL_elem_Slice,GeoGLL_elem_Slice
  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer (kind=MPI_OFFSET_KIND) :: PosFile, PosFileI, PosFileVelR,PosFileGeoR,seek_offset=0
  integer,dimension(:,:),pointer::NbElem_inSlices

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
     write (fnamef,"(a,I4.4,a)") "SField/VelocSlices",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

     do ns=1,Tdomain%nbSlices

        call MPI_ALLGATHER(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize),1,&
        MPI_INTEGER,NbElem_inSlices(ns,:),1,MPI_INTEGER, MPI_COMM_WORLD,code)
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

              select case(int(Tdomain%Slices(ns,1)))
              case (0)
                 allocate(GeoGLL_elem_Slice(0:nglly-1,0:ngllz-1,2))
                 allocate(VelocGLL_elem_Slice(0:nglly-1,0:ngllz-1,0:2))
                 SZG=nglly*ngllz
                 do nx1=0,nglly-1
                    do nx2=0,ngllz-1
                       ipoint = Tdomain%specel(n)%Iglobnum(0,nx1,nx2)

                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (1,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (2,ipoint) 
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(0,nx1,nx2,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(0,nx1,nx2,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(0,nx1,nx2,2)

                    enddo
                 enddo

              case (1)
                 allocate(GeoGLL_elem_Slice(0:ngllx-1,0:ngllz-1,2))
                 allocate(VelocGLL_elem_Slice(0:ngllx-1,0:ngllz-1,0:2))
                 SZG=ngllx*ngllz
                 do nx1=0,ngllx-1
                    do nx2=0,ngllz-1
                       ipoint = Tdomain%specel(n)%Iglobnum(nx1,0,nx2)

                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (0,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (2,ipoint) 
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(nx1,0,nx2,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(nx1,0,nx2,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(nx1,0,nx2,2)

                    enddo
                 enddo

              case (2)

                 allocate(GeoGLL_elem_Slice(0:ngllx-1,0:nglly-1,2))
                 allocate(VelocGLL_elem_Slice(0:ngllx-1,0:nglly-1,0:2))
                 SZG=nglly*ngllx
                 do nx1=0,ngllx-1
                    do nx2=0,nglly-1
                       ipoint = Tdomain%specel(n)%Iglobnum(nx1,nx2,ngllz-1)
                       GeoGLL_elem_Slice(nx1,nx2,1)=Tdomain%GlobCoord (0,ipoint) 
                       GeoGLL_elem_Slice(nx1,nx2,2)=Tdomain%GlobCoord (1,ipoint) 
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,2)
                    enddo
                 enddo
              endselect
              PosFileGeoR=(ne_above+ne_left)*NbOctReal*SZG*2
              PosFileVelR=(ne_above+ne_left)*NbOctReal*SZG*3
              call MPI_FILE_WRITE_AT(descGeo,PosFileGeoR,GeoGLL_elem_Slice(:,:,:),2*SZG,MPI_DOUBLE_PRECISION,statut,code)
              call MPI_FILE_WRITE_AT(desc,PosFileVelR,VelocGLL_elem_Slice(:,:,:),3*SZG,MPI_DOUBLE_PRECISION,statut,code)           
              if (Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize).gt.0) deallocate(GeoGLL_elem_Slice,VelocGLL_elem_Slice)
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
     write (fnamef,"(a,I4.4,a)") "SField/VelocSlices",icount,".out"     
     call MPI_FILE_OPEN(MPI_COMM_WORLD,fnamef,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,desc,code)

     do ns=1,Tdomain%nbSlices

        call MPI_ALLGATHER(Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize),1,&
        MPI_INTEGER,NbElem_inSlices(ns,:),1,MPI_INTEGER, MPI_COMM_WORLD,code)
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

              select case(int(Tdomain%Slices(ns,1)))
              case (0)
                 
                 allocate(VelocGLL_elem_Slice(0:nglly-1,0:ngllz-1,0:2))
                 SZG=nglly*ngllz
                 do nx1=0,nglly-1
                    do nx2=0,ngllz-1
                       
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(0,nx1,nx2,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(0,nx1,nx2,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(0,nx1,nx2,2)
                       !print *,'slices_check1',SZG
                    enddo
                 enddo

              case (1)
                 
                 allocate(VelocGLL_elem_Slice(0:ngllx-1,0:ngllz-1,0:2))
                 SZG=ngllx*ngllz
                 do nx1=0,ngllx-1
                    do nx2=0,ngllz-1
                      
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(nx1,0,nx2,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(nx1,0,nx2,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(nx1,0,nx2,2)

                    enddo
                 enddo

              case (2)

                 allocate(VelocGLL_elem_Slice(0:ngllx-1,0:nglly-1,0:2))
                 SZG=nglly*ngllx
                 do nx1=0,ngllx-1
                    do nx2=0,nglly-1
                       VelocGLL_elem_Slice(nx1,nx2,0)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,0)
                       VelocGLL_elem_Slice(nx1,nx2,1)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,1)
                       VelocGLL_elem_Slice(nx1,nx2,2)=Tdomain%specel(n)%Veloc(nx1,nx2,ngllz-1,2)
                    enddo
                 enddo
              endselect
              
              PosFileVelR=(ne_above+ne_left)*NbOctReal*SZG*3
              call MPI_FILE_WRITE_AT(desc,PosFileVelR,VelocGLL_elem_Slice(:,:,:),3*SZG,MPI_DOUBLE_PRECISION,statut,code)           
              if (Tdomain%I_elem_inSlices(ns,Tdomain%MaxiSliceSize).gt.0) deallocate(VelocGLL_elem_Slice)
           enddo
        endif
     enddo
     call MPI_BARRIER(MPI_COMM_WORLD,code)
     call MPI_FILE_CLOSE(desc,code)
    
  endif

deallocate(NbElem_inSlices)

  return

end subroutine save_slices


