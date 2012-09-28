program Cubit2Spec

! Modified by Q.A. TA on 23/march/07 : L82;L1398-1407
implicit none

type :: processor
integer, dimension(:,:), pointer :: E,V
integer, dimension(:), pointer :: Obj
end type

type :: souvenir
type(processor), dimension(:), pointer :: rank 
end type

type(souvenir), dimension(:), pointer :: memory
real :: x_len, y_len, z_len, dx, dy, dz, R1, R2, R, X, Y, D, rcount
real, dimension (:), allocatable :: xco, yco, zco
real, dimension (:,:), allocatable :: Gcoord
integer :: i,j,k,l, i_count, nel, n_elem, n_elem_xy, n_points, n_pts_xy, nx,ny,nz, edgecut, nparts, &
           proc, n, n_vertices, num, ok, n_faces, n_edges, nf,ne,nv, neighbor, neighbor_face, neighbor_edge,n_pointsold
integer :: idummy, n_blocks, num_mat, n_elem_mat, icount, mapping, store,n_mat
integer, dimension (:), allocatable :: corner, neighbor_corner, orient_corner, counter, &
                                       elmnts, dxadj, dxadjncy, part, vwgt, adjwgt, Material, &
                                       nelem_in_proc, which_nodes, nf_shared, ne_shared, nv_shared
integer, dimension (:,:), allocatable :: elmnts_local, which_elem_in_proc, Ipointer, IpointerN, vertices_shared
integer, dimension (:,:), allocatable :: faces, faces_shared, mapping_faces, mapping_faces_shared
integer, dimension (:,:), allocatable :: edges, edges_shared, mapping_edges, mapping_edges_shared

integer, parameter :: n_dim = 3, n_nods = 8,pi = 3.141592653, &
                      etype = 3, numflag = 0, wgtflag = 0, options = 0
logical :: curve, any_PML, free_surface_pres,okorient,existe
logical, dimension(:), allocatable :: L_Proc
character (len=1) :: yes_or_not
character (len=20) :: meshfilename

! Super Object variables
integer, parameter :: N_MAX_POINTS = 10000
integer :: j1,j2,ok1,ok2,nn,down_mat,up_mat,nnv,nne,nnel,np,nes,nvs,dcount,ucount,n_down,ne_down
integer :: answertype, n_faces_on_interf, n_edges_on_interf, n_vertices_on_interf
integer :: SO_n_faces, SO_n_edges, SO_n_vertices
integer, dimension (:), allocatable :: Vertices_on_interfU, Vertices_on_interfD, Elem_glob2loc,NumProc,nnMax
integer, dimension (:), allocatable :: SO_n_faces_loc, SO_Face_Orient, SO_Edge_Orient, ne_SO_shared, nv_SO_shared, nbFaceSO_Elem,Elem_Ref
integer, dimension (:), allocatable :: corner_d, edges_d, corner_u, edges_u, corner_edge, CountCornerElem
integer, dimension (:,:), allocatable :: Faces_on_interf, Edges_on_interf, CountCorner
integer, dimension (:,:), allocatable :: SO_Elem, SO_Face, SO_Edge, SO_Vertex, SO_Face_NearEdges, SO_Face_NearVertices,SO_Face_NearEdges_Orient 
integer, dimension (:,:), allocatable :: edges_SO_shared, mapping_edges_SO_shared,vertices_SO_shared
integer, dimension (:,:,:), allocatable :: Super_Object_Elem, Super_Object_Orient,CountCorner2Elem
logical :: super_object_present, SO_present_loc
logical, dimension(:), allocatable :: Super_Object, SO_log_Edge, SO_log_Vertex  
character (len=20) :: answer, fnamef
type(processor), dimension(:), pointer :: MemorySO

! Neumann
type(processor), dimension(:), pointer :: MemoryNeu
logical :: Neumann_present,Neu_present_loc
integer :: n_faces_Neumann,n_edges_Neumann,n_vertices_Neumann
integer :: Neu_n_faces, Neu_n_edges, Neu_n_vertices
integer, dimension (:), allocatable :: Vertices_on_Neu,Neumann, Neu_Elem, Neu_n_faces_loc, Neu_Face, Neu_Edge, Neu_Vertex, ne_Neu_shared, nv_Neu_shared 
integer, dimension (:,:), allocatable :: Faces_on_Neumann, Neumann_Glob, Edges_on_Neu
integer, dimension (:,:), allocatable :: Neu_Face_NearEdges, Neu_Face_NearVertices, Neu_Face_NearEdges_Orient 
integer, dimension (:,:), allocatable :: edges_Neu_shared, mapping_edges_Neu_shared, vertices_Neu_shared
logical, dimension(:), allocatable :: Neu_log_Edge, Neu_log_Vertex

! Hill
integer :: nnodes,fact
real :: z0,h0,x0,y0,sigx,sigy,height,dzmax

!Free surface
integer :: n_faces_fsurf,n_vertices_fsurf,diff
integer, dimension (:), allocatable :: vertices_on_fsurf
integer, dimension (:,:), allocatable :: Faces_on_fsurf

!!! No difference is made for now between the vertices and the nodes !!!
!!! Changes will be probably necessary when we introduce 27 control points !!!


open (10,file="cubitmesh",status="old")
read (10,*)    ! *HEADING
read (10,*)                 ! Name of the mesh
read (10,*)                 ! Dimension
read (10,*)    n_blocks     ! Number of materials (blocks)
               n_mat=n_blocks
read (10,*)    n_elem       ! Number of elements
read (10,*)    n_points     ! Number of points
read (10,*)

read (10,*)    ! *NODE
write(*,*) 'Number of Points, Elements, Materials:',n_points, n_elem, n_blocks

allocate (xco(0:n_points-1))
allocate (yco(0:n_points-1))
allocate (zco(0:n_points-1))

do i = 0,n_points-1
   read (10,*)  idummy,xco(i),yco(i),zco(i)
enddo
read (10,*) 


allocate (Ipointer(0:n_nods-1,0:n_elem-1))
allocate (IpointerN(0:n_nods-1,0:n_elem-1))
allocate (elmnts(0:8*n_elem-1))
allocate (Material(0:n_elem-1))

do i=0,n_blocks-1
  read (10,*)    ! *ELEMENT
  read (10,*) num_mat,n_elem_mat
  do n = 0,n_elem_mat-1
    read (10,*)  icount,(Ipointer(j,icount-1),j=0,n_nods-1)
    Material(icount-1) = num_mat - 1
!    call transformation(Ipointer,IpointerN,xco,yco,zco)
  ! Transformations
  if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. & 
       ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
       ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 )  ) then
     okorient = .false.
  else
     okorient = .false.     
  endif

  if ( .not. (okorient) ) then
  if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. & 
       ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
       ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
       ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
       ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) )     then
         IpointerN(0,icount-1) = Ipointer(0,icount-1)
         IpointerN(1,icount-1) = Ipointer(1,icount-1)
         IpointerN(2,icount-1) = Ipointer(5,icount-1)
         IpointerN(3,icount-1) = Ipointer(4,icount-1)
         IpointerN(4,icount-1) = Ipointer(3,icount-1)
         IpointerN(5,icount-1) = Ipointer(2,icount-1)
         IpointerN(6,icount-1) = Ipointer(6,icount-1)
         IpointerN(7,icount-1) = Ipointer(7,icount-1)
         Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
!if (icount==31619) print*,'ori1'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
!if (icount==31619) print*,'ori2'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori3'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori4'
   else if ( ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori5'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori6'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori7'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori8'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori9'
   else if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(1,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(3,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(5,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(7,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori10'
   else if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(1,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(4,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(2,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(7,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori11'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori2'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori13'
   else if ( ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(0,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(2,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(4,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(6,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori14'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(0,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(5,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(3,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(6,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori15'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori16'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori17'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(5,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(3,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori18'
   else if ( ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(2,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(4,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori19'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori20'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori21'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(3,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(5,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori22'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(4,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(2,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori23'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori24'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori25'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori26'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori27'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(5,icount-1)
               IpointerN(3,icount-1) = Ipointer(4,icount-1)
               IpointerN(4,icount-1) = Ipointer(3,icount-1)
               IpointerN(5,icount-1) = Ipointer(2,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori28'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(2,icount-1)
               IpointerN(3,icount-1) = Ipointer(3,icount-1)
               IpointerN(4,icount-1) = Ipointer(4,icount-1)
               IpointerN(5,icount-1) = Ipointer(5,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori29'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori30'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori31'
   else if ( ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(0,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(4,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(2,icount-1)
               IpointerN(6,icount-1) = Ipointer(6,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori32'
   else if ( ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(0,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(7,icount-1)
               IpointerN(3,icount-1) = Ipointer(3,icount-1)
               IpointerN(4,icount-1) = Ipointer(1,icount-1)
               IpointerN(5,icount-1) = Ipointer(5,icount-1)
               IpointerN(6,icount-1) = Ipointer(6,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori33'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori34'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(0,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(7,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(1,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(6,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori35'
   else if ( ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(4,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(3,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(5,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(2,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori36'
   else if ( ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(7,icount-1)
               IpointerN(1,icount-1) = Ipointer(3,icount-1)
               IpointerN(2,icount-1) = Ipointer(0,icount-1)
               IpointerN(3,icount-1) = Ipointer(4,icount-1)
               IpointerN(4,icount-1) = Ipointer(6,icount-1)
               IpointerN(5,icount-1) = Ipointer(2,icount-1)
               IpointerN(6,icount-1) = Ipointer(1,icount-1)
               IpointerN(7,icount-1) = Ipointer(5,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori37'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(3,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(4,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(2,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(5,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori38'
   else if ( ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(4,icount-1)
               IpointerN(1,icount-1) = Ipointer(7,icount-1)
               IpointerN(2,icount-1) = Ipointer(3,icount-1)
               IpointerN(3,icount-1) = Ipointer(0,icount-1)
               IpointerN(4,icount-1) = Ipointer(5,icount-1)
               IpointerN(5,icount-1) = Ipointer(6,icount-1)
               IpointerN(6,icount-1) = Ipointer(2,icount-1)
               IpointerN(7,icount-1) = Ipointer(1,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori40'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori41'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(1,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(6,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(0,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(7,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori42'
   else if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(1,icount-1)
               IpointerN(2,icount-1) = Ipointer(2,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(4,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(7,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori43'
   else if ( ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(1,icount-1)
               IpointerN(2,icount-1) = Ipointer(5,icount-1)
               IpointerN(3,icount-1) = Ipointer(6,icount-1)
               IpointerN(4,icount-1) = Ipointer(3,icount-1)
               IpointerN(5,icount-1) = Ipointer(0,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(7,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori44'
   else if ( ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(2,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(5,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(3,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(4,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori45'
   else if ( ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(6,icount-1)
               IpointerN(1,icount-1) = Ipointer(5,icount-1)
               IpointerN(2,icount-1) = Ipointer(1,icount-1)
               IpointerN(3,icount-1) = Ipointer(2,icount-1)
               IpointerN(4,icount-1) = Ipointer(7,icount-1)
               IpointerN(5,icount-1) = Ipointer(4,icount-1)
               IpointerN(6,icount-1) = Ipointer(0,icount-1)
               IpointerN(7,icount-1) = Ipointer(3,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori46'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(2,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(5,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(3,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(4,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori47'
   else if ( ( xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. & 
             ( xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
             ( yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
             ( zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
               IpointerN(0,icount-1) = Ipointer(5,icount-1)
               IpointerN(1,icount-1) = Ipointer(6,icount-1)
               IpointerN(2,icount-1) = Ipointer(2,icount-1)
               IpointerN(3,icount-1) = Ipointer(1,icount-1)
               IpointerN(4,icount-1) = Ipointer(4,icount-1)
               IpointerN(5,icount-1) = Ipointer(7,icount-1)
               IpointerN(6,icount-1) = Ipointer(3,icount-1)
               IpointerN(7,icount-1) = Ipointer(0,icount-1)
               Ipointer(0:7,icount-1) = IpointerN(0:7,icount-1)
if (icount==31619) print*,'ori48'
   endif
   endif

  enddo
  read (10,*)
enddo


read (10,*) ! Super Object
super_object_present = .false.
read (10,*) answertype
if (answertype == 1) then
   super_object_present = .true.
   read (10,*) answertype  ! (0: Fault, 1:Plane Wave)
   if (answertype == 1) read(10,*) down_mat, up_mat
   read (10,*) n_faces_on_interf
   allocate (Faces_on_interf(0:n_faces_on_interf-1,0:3))
   do i = 0,n_faces_on_interf-1
      read (10,*) icount,(Faces_on_interf(i,j),j=0,3)
   enddo
else
   n_faces_on_interf = 0
endif
print*,'Plane Wave',super_object_present,down_mat, up_mat
read (10,*)
Neumann_present = .false.
read (10,*) ! Neumann   
read (10,*) answertype
if (answertype == 1) then
   Neumann_present = .true.
   read(10,*) down_mat
   read (10,*) n_faces_Neumann  
   allocate (Faces_on_Neumann(0:n_faces_Neumann-1,0:3))
   do i = 0,n_faces_Neumann-1
      read (10,*) icount,(Faces_on_Neumann(i,j),j=0,3)
   enddo
else
   n_faces_Neumann = 0
endif
print*,'Neumann',Neumann_present,down_mat
!read (10,*)
!read (10,*) ! Free Surface
!   read (10,*) n_faces_fsurf
!   allocate (Faces_on_fsurf(0:n_faces_fsurf-1,0:3))
!   do i = 0,n_faces_fsurf-1
!      read (10,*) icount,(Faces_on_fsurf(i,j),j=0,3)
!   enddo
!close(10)

!stop

!! Special treatment for Super Object and Neumann !!
n_edges_on_interf = 0
n_vertices_on_interf = 0
if (super_object_present) then
  allocate (Edges_on_interf(0:4*n_faces_on_interf-1,0:1))
  allocate (Vertices_on_interfU(0:4*n_faces_on_interf-1))
  allocate (Vertices_on_interfD(0:4*n_faces_on_interf-1))
  do i = 0,n_faces_on_interf-1
    Faces_on_interf(i,0:3) = Faces_on_interf(i,0:3)-1
    do j=0,3
      j1 = Faces_on_interf(i,j) 
      if (j==3) then
        j2 = Faces_on_interf(i,0)
      else
        j2 = Faces_on_interf(i,j+1)
      endif
      ok1 = 1
      ok2 = 1
      check1 : do k=0,n_edges_on_interf-1
        if ( (Edges_on_interf(k,0) == j1 .and. Edges_on_interf(k,1) == j2) .or. &
             (Edges_on_interf(k,0) == j2 .and. Edges_on_interf(k,1) == j1) ) then
            ok1 = 0
            exit check1
        endif
      enddo check1
      check2 : do k=0,n_vertices_on_interf-1
        if ( Vertices_on_interfU(k) == j1 ) then
            ok2 = 0
            exit check2
        endif
      enddo check2
      if (ok1==1) then 
        Edges_on_interf(n_edges_on_interf,0) = j1
        Edges_on_interf(n_edges_on_interf,1) = j2
        n_edges_on_interf = n_edges_on_interf + 1
      endif
      if (ok2==1) then 
        Vertices_on_interfU(n_vertices_on_interf) = j1
        n_vertices_on_interf = n_vertices_on_interf + 1
      endif
     enddo 
  enddo
  call sort(Vertices_on_interfU,n_vertices_on_interf)
  write (*,*) "Number of faces, edges, vertices on the SO interface", n_faces_on_interf,n_edges_on_interf,n_vertices_on_interf
endif

if (Neumann_present) then
    n_edges_Neumann = 0
    n_vertices_Neumann = 0
    allocate (Edges_on_Neu(0:4*n_faces_Neumann-1,0:1))
    allocate (Vertices_on_Neu(0:4*n_faces_Neumann-1))
    do i = 0,n_faces_Neumann-1
      Faces_on_Neumann(i,0:3) = Faces_on_Neumann(i,0:3)-1
      do j=0,3
        j1 = Faces_on_Neumann(i,j) 
        if (j==3) then
          j2 = Faces_on_Neumann(i,0)
        else
          j2 = Faces_on_Neumann(i,j+1)
        endif
        ok1 = 1
        ok2 = 1
        check_Neu1 : do k=0,n_edges_Neumann-1
        if ( (Edges_on_Neu(k,0) == j1 .and. Edges_on_Neu(k,1) == j2) .or. &
             (Edges_on_Neu(k,0) == j2 .and. Edges_on_Neu(k,1) == j1) ) then
              ok1 = 0
              exit check_Neu1
        endif
        enddo check_Neu1
        check_Neu2 : do k=0,n_vertices_Neumann-1
          if ( Vertices_on_Neu(k) == j1 ) then
            ok2 = 0
            exit check_Neu2
          endif
        enddo check_Neu2
        if (ok1==1) then 
          Edges_on_Neu(n_edges_Neumann,0) = j1
          Edges_on_Neu(n_edges_Neumann,1) = j2
          n_edges_Neumann = n_edges_Neumann + 1
        endif
        if (ok2==1) then 
          Vertices_on_Neu(n_vertices_Neumann) = j1
          n_vertices_Neumann = n_vertices_Neumann + 1
! Receivers
!          if ( abs(xco(j1)-737000) < 1 ) write(50,*) xco(j1),yco(j1),zco(j1)
!          if ( abs(xco(j1)-730000) < 40 ) write(51,*) xco(j1),yco(j1),zco(j1)          
!
        endif
      enddo 
    enddo
    write (*,*) "Number of faces, edges, vertices on the Neumann boundary", n_faces_Neumann,n_edges_Neumann,n_vertices_Neumann
    deallocate(Edges_on_Neu,Vertices_on_Neu)
endif

! Free surface
!    n_vertices_fsurf = 0
!    allocate (Vertices_on_fsurf(0:4*n_faces_fsurf-1))
!    do i = 0,n_faces_fsurf-1
!      Faces_on_fsurf(i,0:3) = Faces_on_fsurf(i,0:3)-1
!      do j=0,3
!        j1 = Faces_on_fsurf(i,j) 
!        if (j==3) then
!          j2 = Faces_on_fsurf(i,0)
!        else
!          j2 = Faces_on_fsurf(i,j+1)
!        endif
!        ok2 = 1
!        check_fsurf : do k=0,n_vertices_fsurf-1
!          if ( Vertices_on_fsurf(k) == j1 ) then
!            ok2 = 0
!            exit check_fsurf
!          endif
!        enddo check_fsurf
!        if (ok2==1) then 
!          Vertices_on_fsurf(n_vertices_fsurf) = j1
!          n_vertices_fsurf = n_vertices_fsurf + 1
!        endif
!      enddo 
!    enddo
!    write (*,*) "Number of faces, vertices on the Free surface", n_faces_fsurf,n_vertices_fsurf
!    
!open(20,file="freesurf",status="unknown")
!!open(40,file="LPG",status="unknown")
!!open(41,file="SB",status="unknown")
!!open(42,file="EW",status="unknown")
!do n = 0,n_vertices_fsurf-1
!  i = Vertices_on_fsurf(n)
!!  write(21,*) Vertices_on_fsurf(n)
!  write(20,"(3f,a,I)") xco(i),yco(i),zco(i),' ',2
!!  if ( abs(xco(i)-737000)<1 ) write(40,*) xco(i),yco(i),zco(i)
!!  if ( abs(xco(i)-730000)<40 ) write(41,*) xco(i),yco(i),zco(i)
!!  if ( abs(yco(i)-1161800)<50 ) write(42,*) xco(i),yco(i),zco(i)
!enddo
!close(20)
!!close(40)
!!close(41)
!!close(42)
!deallocate(Vertices_on_fsurf,Faces_on_fsurf)
!!end free surface
!stop

! If a super-object is present the number of points should be incremented
allocate(Gcoord(0:n_points+n_vertices_on_interf-1,0:2))
do i=0,n_points-1
  Gcoord(i,0) = xco(i)
  Gcoord(i,1) = yco(i)
  Gcoord(i,2) = zco(i)
enddo
idummy = n_points
do i = 0, n_vertices_on_interf-1 
  Vertices_on_interfD(i) = idummy
  Gcoord(idummy,0:2) = Gcoord(Vertices_on_interfU(i),0:2)
  idummy = idummy + 1
enddo
n_pointsold = n_points
n_points = idummy
if (super_object_present) write (*,*), "Because of doubling, total number of points changed to",n_points
deallocate(xco,yco,zco)


! Changing numbering (-1) 
do n = 0,n_elem-1
!  write(11,*) n
!  write(11,*) (Ipointer(j,n),j=0,n_nods-1)
  do j = 0,n_nods-1
       Ipointer(j,n) = Ipointer(j,n) - 1
  enddo
enddo


! For Metis
do n = 0, n_elem-1
  do i= 0,n_nods-1
    idummy = n*n_nods+i 
    elmnts(idummy) = Ipointer(i,n)
  enddo
enddo
deallocate(Ipointer,IpointerN)


curve = .false.

! Introduction of the number of parts
inquire (file="N_mpi",EXIST=existe)
if ( existe ) then  
    open (15,file="N_mpi", status="old")
    read (15,*) nparts
    close (15)
else
    write (*,*) "File N_mpi not found, introduce the number of parts to partition the mesh"
    read (*,*) nparts
endif


allocate (vwgt(0:n_elem-1))
vwgt = 0
!!! Partitioning the domain by using METIS !!!
allocate (dxadj(0:n_elem))
allocate (dxadjncy(0:8*n_elem-1))
call METIS_MeshToDual (n_elem, n_points, elmnts, etype, numflag, dxadj, dxadjncy)
allocate (adjwgt(0:dxadj(n_elem)-1))
adjwgt = 0
allocate (part(0:n_elem-1))
if (nparts==1) then
   do nel = 0,n_elem-1
      part(nel) = 0
   enddo
else
   call METIS_PartGraphKway (n_elem, dxadj, dxadjncy(0:dxadj(n_elem)-1), vwgt, adjwgt, &
                          wgtflag, numflag, nparts, options, edgecut, part)
endif

!do n=0,n_elem-1
! if ( Material(n)==0 .or. Material(n)==1 ) then
!   part(n) = 0
! else
!   part(n) = 1
! endif
!enddo
!stop
!!! Now there will be two different numberings: !!!
!!! a local one (i.e. relative to a processor) and a global one (i.e. which concerns the whole domain) !!!
!n_points = n_pointsbis
! Super object global informations


if (super_object_present) then                    ! Properties from the global numbering 
  allocate (corner(0:3))                          ! Glob num of 4 vertices of down face 
  allocate (neighbor_corner(0:3))                 ! Glob num of 4 vertices of opposite up face
  allocate (corner_d(0:3))                        ! Internal num (0:7) of 4 vertices of down face
  allocate (edges_d(0:3))                         ! Internal num (0:11) of 4 edges of down face
  allocate (corner_u(0:3))                        ! Internal num (0:7) of 4 vertices of opposite up face
  allocate (edges_u(0:3))                         ! Internal num (0:11) of 4 edges of opposite up face
  allocate (corner_edge(0:1))                     
  allocate (Super_Object(0:n_elem-1))             ! True if an element has a face on the interface
  allocate (Super_Object_Elem(0:n_elem-1,0:3,0:9))    
  allocate (Super_Object_Orient(0:n_elem-1,0:3,0:4))
  allocate (nbFaceSO_Elem(0:n_elem-1))
  Super_Object = .false.
  nbFaceSO_Elem(0:n_elem-1) = 0
  do n = 0, n_elem-1
    if (Material(n) == down_mat) then
      ! Put opposite elements in the same proc and search for opposite faces, edges and vertices
      do nf = 0,5
        select case (nf)  
         case (0)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+1)
          corner(2) = elmnts(8*n+2)
          corner(3) = elmnts(8*n+3)
         case (1)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+1)
          corner(2) = elmnts(8*n+5)
          corner(3) = elmnts(8*n+4)
         case (2)
          corner(0) = elmnts(8*n+1)
          corner(1) = elmnts(8*n+2)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+5)
         case (3)
          corner(0) = elmnts(8*n+3)
          corner(1) = elmnts(8*n+2)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+7)
         case (4)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+3)
          corner(2) = elmnts(8*n+7)
          corner(3) = elmnts(8*n+4)
         case (5)
          corner(0) = elmnts(8*n+4)
          corner(1) = elmnts(8*n+5)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+7)
        end select 
        call caract_face(nf,corner_d,edges_d)  
        search0 : do nn = 0,n_elem-1
           if (Material(nn) == up_mat) then
             search1 : do j = 0,3   ! DOES THE FACE BELONG TO THIS ELEMENT ? 
                 num = corner(j)
                 ok = 0
                 search2 : do k = 0,7
                     if (elmnts(8*nn+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit search2
                     endif
                 enddo search2
                 if (ok==0) exit search1   ! NO, so let's see another element
                 if (j==3) then   ! YES
                   corner_u(0:3) = neighbor_corner(0:3)
                   call sort(neighbor_corner,4)
                   if (neighbor_corner(0)==0) then   ! So which face of the neighbor is it ?
                   if (neighbor_corner(3)==3) neighbor_face = 0
                   if (neighbor_corner(3)==5) neighbor_face = 1
                   if (neighbor_corner(3)==7) neighbor_face = 4
                   else if (neighbor_corner(0)==1) then
                         neighbor_face = 2
                   else if (neighbor_corner(0)==2) then
                         neighbor_face = 3
                   else if (neighbor_corner(0)==4) then
                         neighbor_face = 5
                   else
                         print *,"Coherency Pb between faces and nodes of an element"
                   endif
                   if (part(n) /= part(nn)) part(n) = part (nn)
                   Super_Object(n) = .true.
                   Super_Object(nn) = .true.
                   Super_Object_Elem(nn,nbFaceSO_Elem(nn),0) = n                  ! For an up element: opposite element in down
                   Super_Object_Elem(nn,nbFaceSO_Elem(nn),1) = neighbor_face      ! For an up element: number of its face on the interface
                   Super_Object_Elem(n,nbFaceSO_Elem(n),0) = nn                   ! For a down element: opposite element in up
                   Super_Object_Elem(n,nbFaceSO_Elem(n),1) = nf                   ! For a down element: number of its face on the interface
                   call face_orientation(corner_u,4,neighbor_face,Super_Object_Orient(nn,nbFaceSO_Elem(nn),0))   ! Orientation of the down face in ref to the up
!                   call which_edge(corner_u,edges_u)
                   ! Reference is up face !!
                   call caract_face(neighbor_face,corner_u,edges_u)
                   do l = 0,3
                     do k = 0,7
                       if (elmnts(8*n+k)==elmnts(8*nn+corner_u(l))) then
                         corner_d(l) = k
                       endif
                     enddo
                   enddo 
                   call which_edge(corner_d,edges_d)
                   Super_Object_Elem(nn,nbFaceSO_Elem(nn),2:5) = edges_u(0:3)     ! For an up element: number of its edges on the interface
                   Super_Object_Elem(nn,nbFaceSO_Elem(nn),6:9) = corner_u(0:3)    ! For an up element: number of its vertices on the interface
                   Super_Object_Elem(n,nbFaceSO_Elem(n),2:5) = edges_d(0:3)       ! For a down element: number of its edges on the interface
                   Super_Object_Elem(n,nbFaceSO_Elem(n),6:9) = corner_d(0:3)      ! For a down element: number of its vertices on the interface
                   do i = 0,3
                     neighbor_edge = edges_d(i)
                     corner_edge(0) = corner_d(i)
                     if (i==2) then
                       corner_edge(0) = corner_d(3)
                       corner_edge(1) = corner_d(2)
                     else if (i==3) then
                       corner_edge(0) = corner_d(0)
                       corner_edge(1) = corner_d(3)
                     else
                       corner_edge(1) = corner_d(i+1)
                     endif
                     call edge_orientation(corner_edge,2,neighbor_edge,Super_Object_Orient(nn,nbFaceSO_Elem(nn),1+i)) ! Orientation of the down edges
                   enddo

                   nbFaceSO_Elem(nn) = nbFaceSO_Elem(nn)+1
                 endif
              enddo search1
            endif
         enddo search0         
      enddo

    endif
  enddo
  deallocate (corner,neighbor_corner,corner_edge,corner_d,corner_u,edges_d,edges_u)
  if (nbFaceSO_Elem(nn)>1) print*,'Plusieurs faces PW dans un elem',nn,nbFaceSO_Elem(nn)
endif

if (Neumann_present) then                    ! Properties from the global numbering 
  allocate (Neumann(0:n_elem-1))
  allocate (Neumann_Glob(0:n_elem-1,0:8))
  allocate (corner(0:3))                          ! Glob num of 4 vertices of element face 
  allocate (corner_d(0:3))                        ! Internal num (0:7) of 4 vertices of element face
  allocate (edges_d(0:3))                         ! Internal num (0:11) of 4 edges of element face
  Neumann = .false.

  do n = 0,n_faces_Neumann-1  
     search_neu0 : do nn = 0,n_elem-1
!      if (Material(nn) == down_mat) then
       do nf = 0,5
        select case (nf)  
         case (0)
          corner(0) = elmnts(8*nn)
          corner(1) = elmnts(8*nn+1)
          corner(2) = elmnts(8*nn+2)
          corner(3) = elmnts(8*nn+3)
         case (1)
          corner(0) = elmnts(8*nn)
          corner(1) = elmnts(8*nn+1)
          corner(2) = elmnts(8*nn+5)
          corner(3) = elmnts(8*nn+4)
         case (2)
          corner(0) = elmnts(8*nn+1)
          corner(1) = elmnts(8*nn+2)
          corner(2) = elmnts(8*nn+6)
          corner(3) = elmnts(8*nn+5)
         case (3)
          corner(0) = elmnts(8*nn+3)
          corner(1) = elmnts(8*nn+2)
          corner(2) = elmnts(8*nn+6)
          corner(3) = elmnts(8*nn+7)
         case (4)
          corner(0) = elmnts(8*nn)
          corner(1) = elmnts(8*nn+3)
          corner(2) = elmnts(8*nn+7)
          corner(3) = elmnts(8*nn+4)
         case (5)
          corner(0) = elmnts(8*nn+4)
          corner(1) = elmnts(8*nn+5)
          corner(2) = elmnts(8*nn+6)
          corner(3) = elmnts(8*nn+7)
        end select 
        call caract_face(nf,corner_d,edges_d)

        search_neu1 : do j = 0,3   
            num = Faces_on_Neumann(n,j)
            ok = 0
            search_neu2 : do k = 0,3
                if (corner(k)==num) then
                   ok = 1
                   exit search_neu2
                endif
            enddo search_neu2
            if (ok==0) exit search_neu1 
            if (j==3) then  
              Neumann(nn) = .true.
              Neumann_Glob(nn,0) = nf
              Neumann_Glob(nn,1:4) = edges_d(0:3)
              Neumann_Glob(nn,5:8) = corner_d(0:3)
            endif
        enddo search_neu1 
        if (j==3) exit search_neu0
       enddo
!      endif
     enddo search_neu0
  enddo
  deallocate (corner,corner_d,edges_d)
endif


if (super_object_present) then

  allocate (CountCorner(0:n_elem-1,0:99)) 
  allocate (CountCornerElem(0:99)) 
  allocate (CountCorner2Elem(0:n_elem-1,0:99,0:1))
  CountCorner = 0
  CountCornerElem = 0
  CountCorner2Elem = -1
  allocate (corner(0:3))                          ! Glob num of 4 vertices of down face    
  dcount = 0                 
  do n = 0, n_elem-1
    if ( Material(n) == down_mat .and. (.not. Super_Object(n)) ) then

      ucount = 0
      do nf = 0,5
        select case (nf)  
         case (0)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+1)
          corner(2) = elmnts(8*n+2)
          corner(3) = elmnts(8*n+3)
         case (1)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+1)
          corner(2) = elmnts(8*n+5)
          corner(3) = elmnts(8*n+4)
         case (2)
          corner(0) = elmnts(8*n+1)
          corner(1) = elmnts(8*n+2)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+5)
         case (3)
          corner(0) = elmnts(8*n+3)
          corner(1) = elmnts(8*n+2)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+7)
         case (4)
          corner(0) = elmnts(8*n)
          corner(1) = elmnts(8*n+3)
          corner(2) = elmnts(8*n+7)
          corner(3) = elmnts(8*n+4)
         case (5)
          corner(0) = elmnts(8*n+4)
          corner(1) = elmnts(8*n+5)
          corner(2) = elmnts(8*n+6)
          corner(3) = elmnts(8*n+7)
        end select 

        search10 : do nn = 0,n_elem-1
           if ( Material(nn) == up_mat .and. Super_Object(nn) ) then
             icount = 0
             search11 : do j = 0,3   ! DOES THE FACE BELONG TO THIS ELEMENT ? 
                 num = corner(j)
                 search12 : do k = 0,7
                     if (elmnts(8*nn+k)==num) then
                         icount = icount+1
                         exit search12
                     endif
                 enddo search12
                 if ( j==3 .and. (icount==2 .or. icount==1) )  then
                    if ( icount > CountCorner(dcount,ucount))  then
                      CountCorner(dcount,ucount) = icount
                      CountCorner2Elem(dcount,ucount,0) = n
                      CountCorner2Elem(dcount,ucount,1) = nn 
                    endif
                    ucount = ucount + 1
                 endif
              enddo search11
            endif
         enddo search10         
      enddo
    dcount = dcount+1
    endif
  enddo
  deallocate (corner)
  allocate (nnMax(0:0))
  do n = 0, dcount-1
    CountCornerElem(:) = CountCorner(n,:)
    nnMax = maxloc(CountCornerElem)
    if ( CountCornerElem(nnMax(0)-1) > 0 ) then
      if ( part(CountCorner2Elem(n,nnMax(0)-1,0)) /= part(CountCorner2Elem(n,nnMax(0)-1,1)) )  &
             part(CountCorner2Elem(n,nnMax(0)-1,0)) = part(CountCorner2Elem(n,nnMax(0)-1,1))
    endif
  enddo  
  deallocate(CountCorner,CountCornerElem,CountCorner2Elem)

  do n = 0, n_elem-1
    if ( Material(n)==down_mat ) then  !.and. (.not. Super_Object(n))
      ! Change indices in the down part
      do k = 0,n_nods-1
          do_nn_interior : do  nn = 0, n_vertices_on_interf-1 
             idummy = n*n_nods+k
	     if (elmnts(idummy) == Vertices_on_interfU(nn)) then
	           elmnts(idummy) = Vertices_on_interfD(nn)
!                   if ( part(n) /= NumProc(nn) .and. (.not. Super_Object(n)) ) then
!                      print*,n,part(n),NumProc(nn)
!                      part(n) = NumProc(nn)
!                   endif
		   exit do_nn_interior
	     endif
	  enddo do_nn_interior
      enddo
    endif
  enddo

  ! Change indices in Neumann faces
  if (Neumann_present) then
    do nf=0,n_faces_Neumann-1
      do i=0,3
        do_Neu : do j=0,n_vertices_on_interf-1
          if (Faces_on_Neumann(nf,i)==Vertices_on_interfU(j)) then
            Faces_on_Neumann(nf,i) = Vertices_on_interfD(j)
            exit do_Neu
          endif
        enddo do_Neu
      enddo
    enddo
  endif
endif

! Add vertices in Free surface at the interface
!diff = n_points-n_pointsold
!do i=0,diff-1
!   Vertices_on_fsurf(n_vertices_fsurf+i) = n_points+i
!enddo




!!! Defining the number of elements per processor !!!
allocate (nelem_in_proc(0:nparts-1))
nelem_in_proc = 0
do nel = 0,n_elem-1
    nelem_in_proc(part(nel)) = nelem_in_proc(part(nel)) + 1
enddo

!!! Defining for each processor which elements are inside !!!
!!! "Counter" refers to the local numberings and "nel" to the global numbering !!!
allocate (counter(0:nparts-1))
counter = 0
allocate (which_elem_in_proc(0:nparts-1,0:maxval(nelem_in_proc)-1))
allocate (Elem_glob2loc(0:n_elem-1))                                 ! gives the local number from the global number of an element
allocate (SO_n_faces_loc(0:nparts-1))                                ! number of so faces in each proc
allocate (Neu_n_faces_loc(0:nparts-1))                               ! number of Neumann faces in each proc
SO_n_faces_loc = 0
Neu_n_faces_loc = 0
do nel = 0,n_elem-1
    num = part(nel)
    which_elem_in_proc(num,counter(num)) = nel
    Elem_glob2loc(nel) = counter(num) 
    counter(num) = counter(num) + 1
    if ( super_object_present ) then
      if ( Super_Object(nel) .and. Material(nel) == up_mat )  SO_n_faces_loc(num) = SO_n_faces_loc(num) + nbFaceSO_Elem(nel)
    endif
    if ( Neumann_present ) then
      if ( Neumann(nel) ) Neu_n_faces_loc(num) = Neu_n_faces_loc(num) + 1
    endif
enddo
deallocate (counter)

!print*,'SO_n_faces_loc',SO_n_faces_loc
!print*,'Neu_n_faces_loc',Neu_n_faces_loc


!!! Allocating "memory" !!!
!!! This variable will allow to establish the correspondence between objects shared by different processors !!!
!rcount = 0
allocate (memory(0:n_elem-1))

do nel = 0,n_elem-1
    if (part(nel) /= nparts-1) then
        allocate (memory(nel)%rank(part(nel)+1:nparts-1))
        do proc = part(nel)+1,nparts-1
!            allocate (memory(nel)%rank(proc)%E(0:nelem_in_proc(proc)-1,0:12))
            allocate (memory(nel)%rank(proc)%Obj(0:17))
        enddo
    endif
enddo

if (super_object_present) then
  allocate (MemorySO(0:nparts-1))
  do n = 0,nparts-1
!    allocate (MemorySO(n)%E(0:nparts-1,0:n_points-1))
    allocate (MemorySO(n)%E(0:nparts-1,0:10000))
    MemorySO(n)%E = -1
!    allocate (MemorySO(n)%V(0:nparts-1,0:n_points-1))
    allocate (MemorySO(n)%V(0:nparts-1,0:10000))
    MemorySO(n)%V = -1
  enddo
endif

if (Neumann_present) then
  allocate (MemoryNeu(0:nparts-1))
  do n = 0,nparts-1
!    allocate (MemoryNeu(n)%E(0:nparts-1,0:n_points-1))
    allocate (MemoryNeu(n)%E(0:nparts-1,0:10000))
    MemoryNeu(n)%E = -1
!    allocate (MemoryNeu(n)%V(0:nparts-1,0:n_points-1))
    allocate (MemoryNeu(n)%V(0:nparts-1,0:10000))
    MemoryNeu(n)%V = -1
  enddo
endif


!!! LET'S CONSIDER A PROCESSOR AND CREATE THE FIELDS REQUIRED TO BUILD ITS MESHFILE !!!
call system ("rm -f mesh4spec.????")
meshfilename(1:10) = "mesh4spec."
do proc = 0,nparts-1
 
 if ( SO_n_faces_loc(proc) == 0 ) then 
    SO_present_loc = .false.
 else
    SO_present_loc = .true.
 endif
 if ( Neu_n_faces_loc(proc) == 0 ) then
    Neu_present_loc = .false.
 else
    Neu_present_loc = .true.
 endif

 !!! Defining the nodes which belong to the processor !!!
 !!! These nodes will be sorted according to the global numbering !!!
 allocate (which_nodes(0:8*nelem_in_proc(proc)-1))
 n_vertices = 0
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         ok = 1
         check : do j = 0,n_vertices-1
             if (which_nodes(j)==num) then
                 ok = 0
                 exit check
             endif
         enddo check
         if (ok==1) then
             which_nodes(n_vertices) = num
             n_vertices = n_vertices+1
         endif
     enddo
 enddo
 call sort(which_nodes,n_vertices)


 !!! Associating to each element 8 vertices by using a local numbering !!!
 allocate (elmnts_local(0:nelem_in_proc(proc)-1,0:7))
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         build : do j = 0,n_vertices-1
             if (which_nodes(j)==num) then
                 elmnts_local(n,i) = j
                 exit build
             endif
         enddo build
     enddo
 enddo

 !!! Associating to each element 6 faces by using a local numbering !!!
 !!! Defining the faces shared with another processor !!!
 allocate (faces(0:nelem_in_proc(proc)-1,0:5))
 allocate (mapping_faces(0:nelem_in_proc(proc)-1,0:5))
 allocate (faces_shared(0:nparts-1,0:6*nelem_in_proc(proc)-1))
 allocate (mapping_faces_shared(0:nparts-1,0:6*nelem_in_proc(proc)-1))
 allocate (nf_shared(0:nparts-1))
 mapping = -1
 mapping_faces = -1
 mapping_faces_shared = -1
 nf_shared = 0
 allocate (corner(0:3))
 allocate (neighbor_corner(0:3))
 allocate (orient_corner(0:3))
 n_faces = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     do nf = 0,5   ! nf indicates which face of the element we're considering
         select case (nf)   ! Here we pick up the vertices (in the global numbering) which define the face
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+2)
          corner(3) = elmnts(8*nel+3)
         case (1)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+5)
          corner(3) = elmnts(8*nel+4)
         case (2)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+5)
         case (3)
          corner(0) = elmnts(8*nel+3)
          corner(1) = elmnts(8*nel+2)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         case (4)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
          corner(2) = elmnts(8*nel+7)
          corner(3) = elmnts(8*nel+4)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         end select
         find0 : do i = dxadj(nel), dxadj(nel+1)-1   ! We look at a neighbor of the element
             neighbor = dxadjncy(i)
             find1 : do j = 0,3   ! DOES THE FACE BELONG TO THIS NEIGHBOR ?
                 num = corner(j)
                 ok = 0
                 find2 : do k = 0,7
                     if (elmnts(8*neighbor+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit find2
                     endif
                 enddo find2
                 if (ok==0) exit find1   ! NO, so let's see another neighbor
                 if (j==3) then   ! YES
                     orient_corner(0:3) = neighbor_corner(0:3)
                     call sort(neighbor_corner,4)
                     if (neighbor_corner(0)==0) then   ! So which face of the neighbor is it ?
                        if (neighbor_corner(3)==3) neighbor_face = 0
                        if (neighbor_corner(3)==5) neighbor_face = 1
                        if (neighbor_corner(3)==7) neighbor_face = 4
                     else if (neighbor_corner(0)==1) then
                        neighbor_face = 2
                     else if (neighbor_corner(0)==2) then
                        neighbor_face = 3
                     else if (neighbor_corner(0)==4) then
                        neighbor_face = 5
                     else
                        print *,"Coherency Pb between faces and nodes of an element"
                     endif
                     if (part(neighbor)==proc) then   ! The neighbor is on the same processor than the element
                         if (neighbor>nel) then   ! The neighbor is an element we've never seen
                             faces(n,nf) = n_faces
                             mapping_faces(n,nf) = 0
                             n_faces = n_faces + 1
                         else   ! The neighbor is an element we've ever seen
                             g2l : do i_count = 0,n-1
                                 if (which_elem_in_proc(proc,i_count)==neighbor) then
                                     faces(n,nf) = faces(i_count,neighbor_face)
                                     exit g2l
                                 endif
                             enddo g2l 
                             !! Coherency  
                             call face_orientation(orient_corner,4,neighbor_face,mapping_faces(n,nf))
                          endif   
                     else   ! The neighbor is not on the same processor than the element
                         faces(n,nf) = n_faces
                         mapping_faces(n,nf) = 0  
                         num = part(neighbor)
                         !!! Ensuring the correspondence between the faces shared !!!
                         if (num<proc) then   ! We've ever seen the processor of the neighbor
!                             faces_shared(num,memory(neighbor)%rank(proc)%E(n,12)) = n_faces
                              faces_shared(num,memory(neighbor)%rank(proc)%Obj(neighbor_face)) = n_faces
                              !! Coherency  
!                             call face_orientation(orient_corner,4,neighbor_face,&
!                                                   mapping_faces_shared(num,memory(neighbor)%rank(proc)%E(n,12)))
                              call face_orientation(orient_corner,4,neighbor_face,&
                                                   mapping_faces_shared(num,memory(neighbor)%rank(proc)%Obj(neighbor_face)))

                         else   ! We've never seen the processor of the neighbor
                             faces_shared(num,nf_shared(num)) = n_faces
                             memory(nel)%rank(num)%Obj(nf) = nf_shared(num)
!                             glob2loc : do i_count = 0,nelem_in_proc(num)-1
!                                 if (which_elem_in_proc(num,i_count)==neighbor) then
!                                     memory(nel)%rank(num)%E(i_count,12) = nf_shared(num)
!                                     exit glob2loc
!                                 endif
!                             enddo glob2loc 
                             !! Coherency  
                             call face_orientation(orient_corner,4,neighbor_face,&
                                                   mapping_faces_shared(num,nf_shared(num)))
                         endif
                         nf_shared(num) = nf_shared(num) + 1
                         n_faces = n_faces + 1
                     endif
                     exit find0
                 endif
             enddo find1
         enddo find0
         if (ok==0) then   ! The face is not shared by a neighbor
             faces(n,nf) = n_faces
             mapping_faces(n,nf) = 0
             n_faces = n_faces + 1
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner, orient_corner)

 !!! Associating to each element 12 edges by using a local numbering !!!
 !!! Defining the edges shared with others processors !!!
 allocate (edges(0:nelem_in_proc(proc)-1,0:11))
!
 allocate(Elem_Ref(0:12*nelem_in_proc(proc)-1))
!
 allocate (mapping_edges(0:nelem_in_proc(proc)-1,0:11))
 allocate (L_Proc(0:nparts-1))
 allocate (edges_shared(0:nparts-1,0:12*nelem_in_proc(proc)-1))
 allocate (mapping_edges_shared(0:nparts-1,0:12*nelem_in_proc(proc)-1))
 allocate (ne_shared(0:nparts-1))
 mapping = -1
 mapping_edges = -1
 mapping_edges_shared = -1
 ne_shared = 0
 allocate (corner(0:1))
 allocate (neighbor_corner(0:1))
 allocate (orient_corner(0:1))
 n_edges = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     do ne = 0,11   ! ne indicates which edge of the element we're considering
         select case (ne)   ! Here we pick up the vertices (in the global numbering) which define the edge
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
         case (1)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
         case (2)
          corner(0) = elmnts(8*nel+3)
          corner(1) = elmnts(8*nel+2)
         case (3)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
         case (4)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+5)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
         case (6)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+4)
         case (7)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+6)
         case (8)
          corner(0) = elmnts(8*nel+5)
          corner(1) = elmnts(8*nel+6)
         case (9)
          corner(0) = elmnts(8*nel+7)
          corner(1) = elmnts(8*nel+6)
         case (10)
          corner(0) = elmnts(8*nel+3)
          corner(1) = elmnts(8*nel+7)
         case (11)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+7)
         end select
         findbis0 : do i = 0,n-1   ! HAVE WE EVER SEEN THE EDGE IN ANOTHER ELEMENT OF THE SAME PROCESSOR BEFORE?
             neighbor = which_elem_in_proc(proc,i)
             findbis1 : do j = 0,1
                 num = corner(j)
                 ok = 0
                 findbis2 : do k = 0,7
                     if (elmnts(8*neighbor+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit findbis2
                     endif
                 enddo findbis2
                 if (ok==0) exit findbis1
                 if (j==1) then   ! YES
                     orient_corner(0:1) = neighbor_corner(0:1)
                     call sort(neighbor_corner,2)
                     if (neighbor_corner(0)==0) then
                      if (neighbor_corner(1)==1) neighbor_edge = 0
                      if (neighbor_corner(1)==3) neighbor_edge = 3
                      if (neighbor_corner(1)==4) neighbor_edge = 6
                     else if (neighbor_corner(0)==1) then
                      if (neighbor_corner(1)==2) neighbor_edge = 1
                      if (neighbor_corner(1)==5) neighbor_edge = 4
                     else if (neighbor_corner(0)==2) then
                      if (neighbor_corner(1)==3) neighbor_edge = 2
                      if (neighbor_corner(1)==6) neighbor_edge = 7
                     else if (neighbor_corner(0)==3) then
                      neighbor_edge = 10
                     else if (neighbor_corner(0)==4) then
                      if (neighbor_corner(1)==5) neighbor_edge = 5
                      if (neighbor_corner(1)==7) neighbor_edge = 11
                     else if (neighbor_corner(0)==5) then
                      neighbor_edge = 8
                     else if (neighbor_corner(0)==6) then
                      neighbor_edge = 9
                     endif
                     edges(n,ne) = edges(i,neighbor_edge)

                     !! Coherency
                     call edge_orientation(orient_corner,2,neighbor_edge,mapping_edges(n,ne))
                     exit findbis0
                 endif
             enddo findbis1
         enddo findbis0
         if (ok==0 .or. n==0) then   ! NO
             edges(n,ne) = n_edges
             mapping_edges(n,ne) = 0
!
             Elem_Ref(n_edges) = n
!
             n_edges = n_edges + 1
             L_Proc = .true.
             L_Proc(proc) = .false.
             do i = 0,n_elem-1   ! IS THE EDGE SHARED WITH AN ELEMENT FROM ANOTHER PROCESSOR?
                 if (L_Proc(part(i))) then
                     findter1 : do j = 0,1
                         num = corner(j)
                         ok = 0
                         findter2 : do k = 0,7
                             if (elmnts(8*i+k)==num) then
                                 neighbor_corner(j) = k
                                 ok = 1
                                 exit findter2
                             endif
                         enddo findter2
                         if (ok==0) exit findter1
                         if (j==1) then   ! YES
                             num = part(i)
                             !!! Ensuring the correspondence between the edges shared !!!
                             if (num<proc) then   ! It deals with a processor we've ever seen
                                 orient_corner(0:1) = neighbor_corner(0:1)
                                 call sort(neighbor_corner,2)
                                 if (neighbor_corner(0)==0) then
                                  if (neighbor_corner(1)==1) neighbor_edge = 0
                                  if (neighbor_corner(1)==3) neighbor_edge = 3
                                  if (neighbor_corner(1)==4) neighbor_edge = 6
                                 else if (neighbor_corner(0)==1) then
                                  if (neighbor_corner(1)==2) neighbor_edge = 1
                                  if (neighbor_corner(1)==5) neighbor_edge = 4
                                 else if (neighbor_corner(0)==2) then
                                  if (neighbor_corner(1)==3) neighbor_edge = 2
                                  if (neighbor_corner(1)==6) neighbor_edge = 7
                                 else if (neighbor_corner(0)==3) then
                                  neighbor_edge = 10
                                 else if (neighbor_corner(0)==4) then
                                  if (neighbor_corner(1)==5) neighbor_edge = 5
                                  if (neighbor_corner(1)==7) neighbor_edge = 11
                                 else if (neighbor_corner(0)==5) then
                                  neighbor_edge = 8
                                 else if (neighbor_corner(0)==6) then
                                  neighbor_edge = 9
                                 endif
!                                 edges_shared(num,memory(i)%rank(proc)%E(n,neighbor_edge)) = edges(n,ne)   
                                 edges_shared(num,memory(i)%rank(proc)%Obj(neighbor_edge+6)) = edges(n,ne)
                                 !! Coherency
                                 call edge_orientation(orient_corner,2,neighbor_edge,&
                                                       mapping_edges_shared(num,memory(i)%rank(proc)%Obj(neighbor_edge+6)))     
!                                 call edge_orientation(orient_corner,2,neighbor_edge,&
!                                                       mapping_edges_shared(num,memory(i)%rank(proc)%E(n,neighbor_edge)))                     
                             else   ! It deals with a processor we've never seen
                                 edges_shared(num,ne_shared(num)) = edges(n,ne)
                                 memory(nel)%rank(num)%Obj(ne+6) = ne_shared(num)
!                                 global2local : do i_count = 0,nelem_in_proc(num)-1
!                                     if (which_elem_in_proc(num,i_count)==i) then
!                                         memory(nel)%rank(num)%E(i_count,ne) = ne_shared(num)
!                                         exit global2local
!                                     endif
!                                 enddo global2local
!                                ! Coherency
                                orient_corner(0:1) = neighbor_corner(0:1)
                                 call sort(neighbor_corner,2)
                                 if (neighbor_corner(0)==0) then
                                  if (neighbor_corner(1)==1) neighbor_edge = 0
                                  if (neighbor_corner(1)==3) neighbor_edge = 3
                                  if (neighbor_corner(1)==4) neighbor_edge = 6
                                 else if (neighbor_corner(0)==1) then
                                  if (neighbor_corner(1)==2) neighbor_edge = 1
                                  if (neighbor_corner(1)==5) neighbor_edge = 4
                                 else if (neighbor_corner(0)==2) then
                                  if (neighbor_corner(1)==3) neighbor_edge = 2
                                  if (neighbor_corner(1)==6) neighbor_edge = 7
                                 else if (neighbor_corner(0)==3) then
                                  neighbor_edge = 10
                                 else if (neighbor_corner(0)==4) then
                                  if (neighbor_corner(1)==5) neighbor_edge = 5
                                  if (neighbor_corner(1)==7) neighbor_edge = 11
                                 else if (neighbor_corner(0)==5) then
                                  neighbor_edge = 8
                                 else if (neighbor_corner(0)==6) then
                                  neighbor_edge = 9
                                 endif
                                 call edge_orientation(orient_corner,2,neighbor_edge,&
                                                       mapping_edges_shared(num,ne_shared(num)))     
                             endif
                             ne_shared(num) = ne_shared(num) + 1
                             L_Proc(num) = .false.
                         endif
                     enddo findter1
                 endif
             enddo
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner, orient_corner)


 !!! Defining vertices shared with others processors !!!
 allocate (vertices_shared(0:nparts-1,0:7*nelem_in_proc(proc)-1))
 allocate (nv_shared(0:nparts-1))
 nv_shared = 0
 do n = 0,n_vertices-1
     L_Proc = .true.
     L_Proc(proc) = .false.
     do i = 0,n_elem-1
         if (L_Proc(part(i))) then
             do k = 0,7
                 if (elmnts(8*i+k)==which_nodes(n)) then
                     num = part(i)
                     !!! The local numbering of the vertices follows the global numbering !!!
                     !!! Using memory to ensure the correspondence between the processors is therefore useless !!!
                     vertices_shared(num,nv_shared(num)) = n
                     nv_shared(num) = nv_shared(num) + 1
                     L_Proc(num) = .false.
                 endif
             enddo
         endif
     enddo
 enddo
 deallocate (L_Proc)

!print*,'hello1'
! Define super Object properties
allocate (ne_SO_shared(0:nparts-1))
allocate (nv_SO_shared(0:nparts-1))
allocate (ne_Neu_shared(0:nparts-1))
allocate (nv_Neu_shared(0:nparts-1))
ne_SO_shared = 0
nv_SO_shared = 0
ne_Neu_shared = 0
nv_Neu_shared = 0

if (SO_present_loc) then

 allocate (SO_Elem(0:nelem_in_proc(proc)-1,0:2))
 allocate (SO_Face(0:SO_n_faces_loc(proc)-1,0:1))
 allocate (SO_Edge(0:4*SO_n_faces_loc(proc)-1,0:1))
 allocate (SO_Vertex(0:4*SO_n_faces_loc(proc)-1,0:1))
 allocate (SO_log_Edge(0:n_edges-1))
 allocate (SO_log_Vertex(0:n_vertices-1))
 allocate (SO_Face_NearEdges(0:4*SO_n_faces_loc(proc)-1,0:3))
 allocate (SO_Face_NearEdges_Orient(0:4*SO_n_faces_loc(proc)-1,0:3))
 allocate (SO_Face_NearVertices(0:4*SO_n_faces_loc(proc)-1,0:3))
 allocate (SO_Face_Orient(0:SO_n_faces_loc(proc)-1))
 allocate (SO_Edge_Orient(0:4*SO_n_faces_loc(proc)-1))
 allocate (edges_SO_shared(0:nparts-1,0:4*SO_n_faces_loc(proc)-1))
 allocate (mapping_edges_SO_shared(0:nparts-1,0:4*SO_n_faces_loc(proc)-1))
 allocate (vertices_SO_shared(0:nparts-1,0:4*SO_n_faces_loc(proc)-1))

!! Verification of the orientation compared to the reference edges
 SO_log_Edge = .false.
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
  nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering

  if ( Super_Object(nel) .and. Material(nel) == up_mat ) then

    do k = 0,nbFaceSO_Elem(nel)-1
      do i=0,3
        ne = Super_Object_Elem(nel,k,2+i)
        if ( .not. SO_log_Edge(edges(n,ne)) ) then   ! Edge never seen
          SO_log_Edge(edges(n,ne)) = .true.
          neighbor = Super_Object_Elem(nel,k,0)
          n_down = Elem_glob2loc(neighbor)
          ne_down = Super_Object_Elem(neighbor,0,2+i)

          if ( (Elem_Ref(edges(n,ne)) == n)  .and. (.not.(Elem_Ref(edges(n_down,ne_down)) == n_down)) ) then
             if ( mapping_edges(n_down,ne_down) == 1 ) then
                 if ( Super_Object_Orient(nel,k,1+i) == 0 ) then
                    Super_Object_Orient(nel,k,1+i) = 1
                 else
                    Super_Object_Orient(nel,k,1+i) = 0                
                 endif
             endif
          endif
          if ( (Elem_Ref(edges(n_down,ne_down)) == n_down)  .and. (.not.(Elem_Ref(edges(n,ne)) == n)) ) then
             if ( mapping_edges(n,ne) == 1 ) then
                 if ( Super_Object_Orient(nel,k,1+i) == 0 ) then
                    Super_Object_Orient(nel,k,1+i) = 1
                 else
                    Super_Object_Orient(nel,k,1+i) = 0                
                 endif
             endif
          endif
          if ( (.not.(Elem_Ref(edges(n,ne)) == n))  .and. (.not.(Elem_Ref(edges(n_down,ne_down)) == n_down)) ) then
             if ( .not. ( mapping_edges(n,ne) == mapping_edges(n_down,ne_down) ) ) then
                 if ( Super_Object_Orient(nel,k,1+i) == 0 ) then
                    Super_Object_Orient(nel,k,1+i) = 1
                 else
                    Super_Object_Orient(nel,k,1+i) = 0                
                 endif
             endif
          endif
!if ( proc==1 .and. n==348 .and. n_down==8 ) print*,'Test Apres',Super_Object_Orient(nel,k,1+i)
        endif
      enddo
    enddo
  endif
 enddo

 SO_Face_NearEdges = -1
 SO_Face_NearVertices = -1
 SO_Face_NearEdges_Orient = -1
 SO_n_faces = 0
 SO_n_edges = 0
 SO_n_vertices = 0
 SO_Elem = -1
 SO_log_Edge = .false.
 SO_log_Vertex = .false.

 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
  nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering

  if ( Super_Object(nel) .and. Material(nel) == up_mat ) then
    do k = 0,nbFaceSO_Elem(nel)-1
      SO_Elem(n,k) = SO_n_faces
      neighbor = Super_Object_Elem(nel,k,0)
      ! Faces
      nf = Super_Object_Elem(nel,k,1)
      SO_Face(SO_n_faces,0) = faces(n,nf)
      nf = Super_Object_Elem(neighbor,0,1)
      SO_Face(SO_n_faces,1) = faces(Elem_glob2loc(neighbor),nf)
      SO_Face_Orient(SO_n_faces) = Super_Object_Orient(nel,k,0)
      !Edges
      do i=0,3
        ne = Super_Object_Elem(nel,k,2+i)
!if (proc==1) print*,SO_log_Edge(1413),n,i
        if ( .not. SO_log_Edge(edges(n,ne)) ) then   ! Edge never seen
          SO_Edge(SO_n_edges,0) = edges(n,ne)
          SO_log_Edge(SO_Edge(SO_n_edges,0)) = .true.
          SO_Face_NearEdges_Orient(SO_n_faces,i) = mapping_edges(n,ne)
          SO_Edge_Orient(SO_n_edges) = Super_Object_Orient(nel,k,1+i)
!if ( proc==1 .and. edges(n,ne)== 1413) print*,'pas deja compte',SO_n_edges,nel,k,i,proc,edges(n,ne),n
          ne = Super_Object_Elem(neighbor,0,2+i)
          SO_Edge(SO_n_edges,1) = edges(Elem_glob2loc(neighbor),ne)
          SO_log_Edge(SO_Edge(SO_n_edges,1)) = .true.
          SO_Face_NearEdges(SO_n_faces,i) = SO_n_edges
!if (proc==2 .and. n==430) print*,SO_n_edges,SO_Edge(SO_n_edges,0),ne,mapping_edges(n,ne),nel,k,i,Super_Object_Elem(1403,0,3)
!if (proc==0) then
!if ( SO_Edge(SO_n_edges,0) == 8230) write(93,*) 'SO_Edge8230',SO_n_edges,n,ne,nel,i
!if ( SO_Edge(SO_n_edges,1) == 27) write(93,*) 'SO_Edge27',SO_n_edges,n,ne,nel,i,neighbor,Elem_glob2loc(neighbor)
!endif
          do np = 0,nparts-1
            do nes = 0,ne_shared(np)-1
              if (edges_shared(np,nes) == SO_Edge(SO_n_edges,0)) then 

              if (  MemorySO(np)%E(proc,nes) < 0 ) then  ! Edge never seen by the other proc np
                 MemorySO(proc)%E(np,nes) = ne_SO_shared(np)
                 edges_SO_shared(np,MemorySO(proc)%E(np,nes)) =  SO_n_edges
                 mapping_edges_SO_shared(np,MemorySO(proc)%E(np,nes)) = mapping_edges_shared(np,nes)
                 ne_SO_shared(np) = ne_SO_shared(np)+1
              else  ! Si oui, il faut garder le meme ordre d'echange
                 edges_SO_shared(np,MemorySO(np)%E(proc,nes)) =  SO_n_edges
                 mapping_edges_SO_shared(np,MemorySO(np)%E(proc,nes)) = mapping_edges_shared(np,nes)
                 ne_SO_shared(np) = ne_SO_shared(np)+1
              endif
!                edges_SO_shared(np,ne_SO_shared(np)) =  SO_n_edges
!                mapping_edges_SO_shared(np,ne_SO_shared(np)) = mapping_edges_shared(np,nes)
!                ne_SO_shared(np) = ne_SO_shared(np)+1
              endif
            enddo
          enddo
          SO_n_edges = SO_n_edges + 1
         else                                         ! Edge already counted
          do nn = 0,n-1
             nnel = which_elem_in_proc(proc,nn)
             do l = 0,nbFaceSO_Elem(nnel)-1
               if ( SO_Elem(nn,l) > -1 ) then
                 do j = 0,3
                   nne = Super_Object_Elem(nnel,l,2+j)
                   if ( edges(nn,nne) == edges(n,ne) ) then 
                     SO_Face_NearEdges(SO_n_faces,i) = SO_Face_NearEdges(SO_Elem(nn,l),j) 
                     SO_Face_NearEdges_Orient(SO_n_faces,i) =  mapping_edges(n,ne)  ! SO_Face_NearEdges_Orient(SO_Elem(nn,l),j)
!if ( proc==1 .and. edges(n,ne)== 1413) print*,'deja compte',SO_Face_NearEdges(SO_Elem(nn,l),j),edges(n,ne),nn,j,mapping_edges(n,ne) 
                     endif
                 enddo
               endif 
             enddo 
          enddo
        endif
      enddo

      !Vertices
      do i=0,3
        nv = Super_Object_Elem(nel,k,6+i)
        if ( .not. SO_log_Vertex(elmnts_local(n,nv)) ) then   ! Vertex never seen
          SO_Vertex(SO_n_vertices,0) = elmnts_local(n,nv)
          SO_log_Vertex(SO_Vertex(SO_n_vertices,0)) = .true.
          nv = Super_Object_Elem(neighbor,0,6+i)
          SO_Vertex(SO_n_vertices,1) = elmnts_local(Elem_glob2loc(neighbor),nv)
          SO_log_Vertex(SO_Vertex(SO_n_vertices,1)) = .true.
          SO_Face_NearVertices(SO_n_faces,i) = SO_n_vertices
          do np = 0,nparts-1
            do nvs = 0,nv_shared(np)-1
              if (vertices_shared(np,nvs) == SO_Vertex(SO_n_vertices,0)) then 

              if (  MemorySO(np)%V(proc,nvs) < 0 ) then
                 MemorySO(proc)%V(np,nvs) = nv_SO_shared(np)
                 vertices_SO_shared(np,MemorySO(proc)%V(np,nvs)) =  SO_n_vertices
                 nv_SO_shared(np) = nv_SO_shared(np)+1
              else
                 vertices_SO_shared(np,MemorySO(np)%V(proc,nvs)) =  SO_n_vertices
                 nv_SO_shared(np) = nv_SO_shared(np)+1
              endif

!                vertices_SO_shared(np,nv_SO_shared(np)) =  SO_n_vertices
!                nv_SO_shared(np) = nv_SO_shared(np)+1
              endif
            enddo
          enddo
          SO_n_vertices = SO_n_vertices + 1
        else                                              ! Vertex already counted
          do nn = 0,n
             nnel = which_elem_in_proc(proc,nn)
               do l = 0,nbFaceSO_Elem(nnel)-1
                if ( SO_Elem(nn,l) > -1 ) then 
                 do j = 0,3
                   nnv = Super_Object_Elem(nnel,l,6+j)
                   if ( elmnts_local(nn,nnv) == elmnts_local(n,nv) ) then 
                     SO_Face_NearVertices(SO_n_faces,i) = SO_Face_NearVertices(SO_Elem(nn,l),j)   
                  endif
                 enddo
               endif
             enddo
          enddo
   
        endif
      enddo
      SO_n_faces = SO_n_faces + 1
    enddo
  endif
 enddo
endif
!print*,'hello2'
! Neumann
if (Neu_present_loc) then
 allocate (Neu_Elem(0:nelem_in_proc(proc)-1))
 allocate (Neu_Face(0:Neu_n_faces_loc(proc)-1))
 allocate (Neu_Edge(0:4*Neu_n_faces_loc(proc)-1))
 allocate (Neu_Vertex(0:4*Neu_n_faces_loc(proc)-1))
 allocate (Neu_log_Edge(0:n_edges-1))
 allocate (Neu_log_Vertex(0:n_vertices-1))
 allocate (Neu_Face_NearEdges(0:4*Neu_n_faces_loc(proc)-1,0:3))
 allocate (Neu_Face_NearEdges_Orient(0:4*Neu_n_faces_loc(proc)-1,0:3))
 allocate (Neu_Face_NearVertices(0:4*Neu_n_faces_loc(proc)-1,0:3))
 allocate (edges_Neu_shared(0:nparts-1,0:4*Neu_n_faces_loc(proc)-1))
 allocate (mapping_edges_Neu_shared(0:nparts-1,0:4*Neu_n_faces_loc(proc)-1))
 allocate (vertices_Neu_shared(0:nparts-1,0:4*Neu_n_faces_loc(proc)-1))
 Neu_Face_NearEdges = -1
 Neu_Face_NearVertices = -1
 Neu_Face_NearEdges_Orient = -1
 Neu_n_faces = 0
 Neu_n_edges = 0
 Neu_n_vertices = 0
 Neu_Elem = -1
 Neu_log_Edge = .false.
 Neu_log_Vertex = .false.

 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
  nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering

  if ( Neumann(nel) ) then
    Neu_Elem(n) = Neu_n_faces
    ! Faces
    nf = Neumann_Glob(nel,0)  !num (0:5) de la face dans l'element global nel
    Neu_Face(Neu_n_faces) = faces(n,nf) !num local au proc 

    !Edges
    do i=0,3
      ne = Neumann_Glob(nel,1+i)
      if ( .not. Neu_log_Edge(edges(n,ne)) ) then   ! Edge never seen
        Neu_Edge(Neu_n_edges) = edges(n,ne)
        Neu_log_Edge(Neu_Edge(Neu_n_edges)) = .true.
        Neu_Face_NearEdges(Neu_n_faces,i) = Neu_n_edges
        Neu_Face_NearEdges_Orient(Neu_n_faces,i) = mapping_edges(n,ne)

        do np = 0,nparts-1
          do nes = 0,ne_shared(np)-1
            if (edges_shared(np,nes) == Neu_Edge(Neu_n_edges)) then 

              if (  MemoryNeu(np)%E(proc,nes) < 0 ) then
                 MemoryNeu(proc)%E(np,nes) = ne_Neu_shared(np)
                 edges_Neu_shared(np,MemoryNeu(proc)%E(np,nes)) =  Neu_n_edges
                 mapping_edges_Neu_shared(np,MemoryNeu(proc)%E(np,nes)) = mapping_edges_shared(np,nes)
                 ne_Neu_shared(np) = ne_Neu_shared(np)+1
              else
                 edges_Neu_shared(np,MemoryNeu(np)%E(proc,nes)) =  Neu_n_edges
                 mapping_edges_Neu_shared(np,MemoryNeu(np)%E(proc,nes)) = mapping_edges_shared(np,nes)
                 ne_Neu_shared(np) = ne_Neu_shared(np)+1
              endif
  
!              edges_Neu_shared(np,ne_Neu_shared(np)) =  Neu_n_edges
!              mapping_edges_Neu_shared(np,ne_Neu_shared(np)) = mapping_edges_shared(np,nes)
!              ne_Neu_shared(np) = ne_Neu_shared(np)+1
            endif
          enddo
        enddo
        Neu_n_edges = Neu_n_edges + 1
       else                                         ! Edge already counted
        do nn = 0,n-1
           nnel = which_elem_in_proc(proc,nn)
           if ( Neu_Elem(nn) > -1 ) then 
             do j = 0,3
               nne = Neumann_Glob(nnel,1+j)
               if ( edges(nn,nne) == edges(n,ne) ) then 
                  Neu_Face_NearEdges(Neu_n_faces,i) = Neu_Face_NearEdges(Neu_Elem(nn),j) 
                  Neu_Face_NearEdges_Orient(Neu_n_faces,i) = mapping_edges(n,ne) ! Neu_Face_NearEdges_Orient(Neu_Elem(nn),j)
               endif
             enddo
           endif 
        enddo
      endif
     enddo

    !Vertices
    do i=0,3
      nv = Neumann_Glob(nel,5+i)
      if ( .not. Neu_log_Vertex(elmnts_local(n,nv)) ) then   ! Vertex never seen
        Neu_Vertex(Neu_n_vertices) = elmnts_local(n,nv)
        Neu_log_Vertex(Neu_Vertex(Neu_n_vertices)) = .true.
        Neu_Face_NearVertices(Neu_n_faces,i) = Neu_n_vertices
        do np = 0,nparts-1
          do nvs = 0,nv_shared(np)-1
            if (vertices_shared(np,nvs) == Neu_Vertex(Neu_n_vertices)) then

              if (  MemoryNeu(np)%V(proc,nvs) < 0 ) then
                 MemoryNeu(proc)%V(np,nvs) = nv_Neu_shared(np)
                 vertices_Neu_shared(np,MemoryNeu(proc)%V(np,nvs)) =  Neu_n_vertices
                 nv_Neu_shared(np) = nv_Neu_shared(np)+1
              else
                 vertices_Neu_shared(np,MemoryNeu(np)%V(proc,nvs)) =  Neu_n_vertices
                 nv_Neu_shared(np) = nv_Neu_shared(np)+1
              endif
            endif
          enddo
        enddo
        Neu_n_vertices = Neu_n_vertices + 1
      else                                              ! Vertex already counted
        do nn = 0,n-1
           nnel = which_elem_in_proc(proc,nn)
            if ( Neu_Elem(nn) > -1 ) then 
             do j = 0,3
               nnv = Neumann_Glob(nnel,5+j)
               if ( elmnts_local(nn,nnv) == elmnts_local(n,nv) ) then 
                  Neu_Face_NearVertices(Neu_n_faces,i) = Neu_Face_NearVertices(Neu_Elem(nn),j)   
               endif
             enddo
           endif 
        enddo

   
      endif
     enddo
    Neu_n_faces = Neu_n_faces + 1

  endif
 enddo
endif


 !!! Writing the meshfile !!!
 write (meshfilename(11:14), '(i4.4)') proc
 open (11, file = trim(meshfilename))
 write (11,"(1i6,a)") n_dim
 write (11,*)
 write (11,"(1i6,a)") n_vertices
 write (11,*) curve
 do n = 0,n_vertices-1
 i_count = which_nodes(n)
 !    write (11,*) xco(i_count),yco(i_count),zco(i_count)
    write (11,*) Gcoord(i_count,0:2)
 enddo
 write (11,*)
 write (11,"(1i6,a)") nelem_in_proc(proc)
 write (11,*)
 write (11,"(1i6,a)") n_mat
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     write (11,"(1i6)") Material(nel)
 enddo
 write (11,*)
 write (11,"(1i6,a)") n_nods
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(8i6)") (elmnts_local(n,i),i=0,n_nods-1)
 enddo
 write (11,*)
 write (11,"(1i6,a)") n_faces
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(6i6)") (faces(n,i),i=0,5)
     write (11,"(6i6)") (mapping_faces(n,i),i=0,5)
 enddo
 write (11,*)
 write (11,"(1i6,a)") n_edges
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(12i6)") (edges(n,i),i=0,11)
     write (11,"(12i6)") (mapping_edges(n,i),i=0,11)
 enddo
 write (11,*)
 write (11,"(1i6,a)") n_vertices
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(8i6)") (elmnts_local(n,i),i=0,7)
 enddo
 if ( .not. SO_present_loc) then
   write (11,*)
   write (11,*) "Super Object"
   write (11,*)  .false.
 else
   write (11,*)
   write (11,*) "Super Object"
   write (11,*)  .true.
   write (11,*) "Faces in Super Object"
   write (11,"(1i6,a)") SO_n_faces_loc(proc)
   write (11,*) "4 Edges for each face"
   do n = 0,SO_n_faces-1
     write (11,"(4i6)") (SO_Face_NearEdges(n,i),i=0,3)
     write (11,"(4i6)") (SO_Face_NearEdges_Orient(n,i),i=0,3)
   enddo
   write (11,*) "4 Vertices for each face"
   do n = 0,SO_n_faces-1
     write (11,"(4i6)") (SO_Face_NearVertices(n,i),i=0,3)
   enddo
   write (11,*) "Glob number of up and down faces and orientation of down compared to up"
   write (11,*) SO_n_faces
   do n = 0,SO_n_faces-1
     write (11,"(3i6)") SO_Face(n,0),SO_Face(n,1),SO_Face_Orient(n)
   enddo
   write (11,*) "Glob number of up and down edges and orientation of down compared to up"
   write (11,*) SO_n_edges
   do n = 0,SO_n_edges-1
     write (11,"(3i6)") SO_Edge(n,0),SO_Edge(n,1),SO_Edge_Orient(n)
   enddo
   write (11,*) "Glob number of up and down vertices"
   write (11,*) SO_n_vertices
   do n = 0,SO_n_vertices-1
     write (11,"(2i6)") SO_Vertex(n,0),SO_Vertex(n,1)
   enddo
 endif
 if ( .not. Neu_present_loc) then
   write (11,*)
   write (11,*) "Neumann"
   write (11,*)  .false.
 else
   write (11,*)
   write (11,*) "Neumann"
   write (11,*)  .true.
   write (11,*) "Faces in Neumann"
   write (11,"(1i6,a)") Neu_n_faces_loc(proc)
   write (11,*) "4 Edges for each face"
   do n = 0,Neu_n_faces-1
     write (11,"(4i6)") (Neu_Face_NearEdges(n,i),i=0,3)
     write (11,"(4i6)") (Neu_Face_NearEdges_Orient(n,i),i=0,3)
   enddo
   write (11,*) "4 Vertices for each face"
   do n = 0,Neu_n_faces-1
     write (11,"(4i6)") (Neu_Face_NearVertices(n,i),i=0,3)
   enddo
   write (11,*) "Glob number of faces"
   write (11,*) Neu_n_faces
   do n = 0,Neu_n_faces-1
     write (11,"(1i6)") Neu_Face(n)
   enddo
   write (11,*) "Glob number of edges"
   write (11,*) Neu_n_edges
   do n = 0,Neu_n_edges-1
     write (11,"(1i6)") Neu_Edge(n)
   enddo
   write (11,*) "Glob number of vertices"
   write (11,*) Neu_n_vertices
   do n = 0,Neu_n_vertices-1
     write (11,"(1i6)") Neu_Vertex(n)
   enddo
 endif
 write (11,*)
 write (11,"(1i6,a)") nparts
 do n = 0,nparts-1
     write (11,"(7i6)") nf_shared(n),ne_shared(n),nv_shared(n),ne_SO_shared(n),nv_SO_shared(n),ne_Neu_shared(n),nv_Neu_shared(n)
     do nf = 0,nf_shared(n)-1
         write (11,"(2i6)") faces_shared(n,nf),mapping_faces_shared(n,nf)
     enddo
     do ne = 0,ne_shared(n)-1
         write (11,"(2i6)") edges_shared(n,ne),mapping_edges_shared(n,ne)
     enddo
     do nv = 0,nv_shared(n)-1
         write (11,"(1i6)") vertices_shared(n,nv)
     enddo
     do ne = 0,ne_SO_shared(n)-1
         write (11,"(2i6)") edges_SO_shared(n,ne),mapping_edges_SO_shared(n,ne)
     enddo
     do nv = 0,nv_SO_shared(n)-1
         write (11,"(1i6)") vertices_SO_shared(n,nv)
     enddo
     do ne = 0,ne_Neu_shared(n)-1
         write (11,"(2i6)") edges_Neu_shared(n,ne),mapping_edges_Neu_shared(n,ne)
     enddo
     do nv = 0,nv_Neu_shared(n)-1
         write (11,"(1i6)") vertices_Neu_shared(n,nv)
     enddo
 enddo
 write (11,*)
 write (11,*) 'Free Surf'
 write (11,*) 'F'
 close (11)

 deallocate (which_nodes, elmnts_local, faces, edges, Elem_Ref)
 deallocate (nf_shared, ne_shared, nv_shared, faces_shared, edges_shared, vertices_shared,ne_Neu_shared,nv_Neu_shared,ne_SO_shared,nv_SO_shared)
 deallocate (mapping_faces, mapping_faces_shared, mapping_edges, mapping_edges_shared)

 if (SO_present_loc) then
   deallocate (SO_Elem,SO_Face,SO_Edge,SO_Vertex,SO_log_Edge,SO_log_Vertex)
   deallocate (SO_Face_NearEdges,SO_Face_NearVertices,SO_Face_Orient,SO_Edge_Orient,SO_Face_NearEdges_Orient)
   deallocate (edges_SO_shared,mapping_edges_SO_shared,vertices_SO_shared)
 endif

 if (Neu_present_loc) then
     deallocate (Neu_Elem,Neu_Face,Neu_Edge,Neu_Vertex,Neu_log_Edge,Neu_log_Vertex)
     deallocate (Neu_Face_NearEdges,Neu_Face_NearEdges_Orient,Neu_Face_NearVertices)
     deallocate (edges_Neu_shared,mapping_edges_Neu_shared,vertices_Neu_shared)
 endif

enddo

!do nel = 0,n_elem-1
!    if (part(nel) /= nparts-1) then
!        do proc = part(nel)+1,nparts-1
!            deallocate (memory(nel)%rank(proc)%E)
!        enddo
!        deallocate (memory(nel)%rank)
!    endif
!enddo

deallocate (elmnts, adjwgt, dxadj, dxadjncy, vwgt, part)
deallocate (Material, nelem_in_proc, which_elem_in_proc, Elem_glob2loc,SO_n_faces_loc,Neu_n_faces_loc)

if (super_object_present) then
   deallocate (Faces_on_interf,Edges_on_interf,Vertices_on_interfU,Vertices_on_interfD)       
   deallocate (Super_Object,Super_Object_Elem,Super_Object_Orient,nbFaceSO_Elem)
endif

if (Neumann_present) then
    deallocate (Faces_on_Neumann,Neumann,Neumann_Glob)
endif

end program Cubit2Spec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sort(vector, length)

implicit none

integer, intent(IN) :: length
integer, dimension(0:length-1), intent(INOUT) :: vector
integer :: n, ok, i, tmp

n = length-1
ok = 1
do while (ok==1)
    ok = 0
    do i = 0,n-1
        if (vector(i)>vector(i+1)) then
            tmp = vector(i+1)
            vector(i+1) = vector(i)
            vector(i) = tmp
            ok = 1
        endif
    enddo
    n = n-1
enddo

return
end subroutine sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine face_orientation(corner,length,nf,orient)

implicit none

integer, intent(IN) :: length, nf
integer, dimension(0:length-1), intent(IN) :: corner
integer, intent(INOUT) :: orient

integer :: n1,n2,n3

n1 = corner(0)
n2 = corner(1)
n3 = corner(2)
if (nf==0) then
  if ( n1==0 .and. n2==1 .and. n3==2 ) orient = 0
  if ( n1==1 .and. n2==2 .and. n3==3 ) orient = 6
  if ( n1==2 .and. n2==3 .and. n3==0 ) orient = 3
  if ( n1==3 .and. n2==0 .and. n3==1 ) orient = 5
  if ( n1==0 .and. n2==3 .and. n3==2 ) orient = 4
  if ( n1==3 .and. n2==2 .and. n3==1 ) orient = 2
  if ( n1==2 .and. n2==1 .and. n3==0 ) orient = 7
  if ( n1==1 .and. n2==0 .and. n3==3 ) orient = 1
else if (nf==1) then
  if ( n1==0 .and. n2==1 .and. n3==5 ) orient = 0
  if ( n1==1 .and. n2==5 .and. n3==4 ) orient = 6
  if ( n1==5 .and. n2==4 .and. n3==0 ) orient = 3
  if ( n1==4 .and. n2==0 .and. n3==1 ) orient = 5
  if ( n1==0 .and. n2==4 .and. n3==5 ) orient = 4
  if ( n1==4 .and. n2==5 .and. n3==1 ) orient = 2
  if ( n1==5 .and. n2==1 .and. n3==0 ) orient = 7
  if ( n1==1 .and. n2==0 .and. n3==4 ) orient = 1
else if (nf==2) then
  if ( n1==1 .and. n2==2 .and. n3==6 ) orient = 0
  if ( n1==2 .and. n2==6 .and. n3==5 ) orient = 6
  if ( n1==6 .and. n2==5 .and. n3==1 ) orient = 3
  if ( n1==5 .and. n2==1 .and. n3==2 ) orient = 5
  if ( n1==1 .and. n2==5 .and. n3==6 ) orient = 4
  if ( n1==5 .and. n2==6 .and. n3==2 ) orient = 2
  if ( n1==6 .and. n2==2 .and. n3==1 ) orient = 7
  if ( n1==2 .and. n2==1 .and. n3==5 ) orient = 1
else if (nf==3) then
  if ( n1==3 .and. n2==2 .and. n3==6 ) orient = 0
  if ( n1==2 .and. n2==6 .and. n3==7 ) orient = 6
  if ( n1==6 .and. n2==7 .and. n3==3 ) orient = 3
  if ( n1==7 .and. n2==3 .and. n3==2 ) orient = 5
  if ( n1==3 .and. n2==7 .and. n3==6 ) orient = 4
  if ( n1==7 .and. n2==6 .and. n3==2 ) orient = 2
  if ( n1==6 .and. n2==2 .and. n3==3 ) orient = 7
  if ( n1==2 .and. n2==3 .and. n3==7 ) orient = 1
else if (nf==4) then
  if ( n1==0 .and. n2==3 .and. n3==7 ) orient = 0
  if ( n1==3 .and. n2==7 .and. n3==4 ) orient = 6
  if ( n1==7 .and. n2==4 .and. n3==0 ) orient = 3
  if ( n1==4 .and. n2==0 .and. n3==3 ) orient = 5
  if ( n1==0 .and. n2==4 .and. n3==7 ) orient = 4
  if ( n1==4 .and. n2==7 .and. n3==3 ) orient = 2
  if ( n1==7 .and. n2==3 .and. n3==0 ) orient = 7
  if ( n1==3 .and. n2==0 .and. n3==4 ) orient = 1
else if (nf==5) then
  if ( n1==4 .and. n2==5 .and. n3==6 ) orient = 0
  if ( n1==5 .and. n2==6 .and. n3==7 ) orient = 6
  if ( n1==6 .and. n2==7 .and. n3==4 ) orient = 3
  if ( n1==7 .and. n2==4 .and. n3==5 ) orient = 5
  if ( n1==4 .and. n2==7 .and. n3==6 ) orient = 4
  if ( n1==7 .and. n2==6 .and. n3==5 ) orient = 2
  if ( n1==6 .and. n2==5 .and. n3==4 ) orient = 7
  if ( n1==5 .and. n2==4 .and. n3==7 ) orient = 1
endif


end subroutine face_orientation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine edge_orientation(corner,length,ne,orient)

implicit none

integer, intent(IN) :: length, ne
integer, dimension(0:length-1), intent(IN) :: corner
integer, intent(INOUT) :: orient

orient = 0
if (ne==2 .or. ne==9) then
  if (corner(0) < corner(1)) orient = 1
else 
  if (corner(0) > corner(1)) orient = 1
endif


end subroutine edge_orientation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine caract_face(nf,corner_d,edge_d)

implicit none

integer, intent(IN) :: nf
integer, dimension(0:3), intent(INOUT) :: corner_d,edge_d


   select case (nf)  
       case (0)
          corner_d(0) = 0
          corner_d(1) = 1
          corner_d(2) = 2
          corner_d(3) = 3
          edge_d(0) = 0
          edge_d(1) = 1
          edge_d(2) = 2
          edge_d(3) = 3
       case (1)
          corner_d(0) = 0
          corner_d(1) = 1
          corner_d(2) = 5
          corner_d(3) = 4
          edge_d(0) = 0
          edge_d(1) = 4
          edge_d(2) = 5
          edge_d(3) = 6
       case (2)
          corner_d(0) = 1
          corner_d(1) = 2
          corner_d(2) = 6
          corner_d(3) = 5
          edge_d(0) = 1
          edge_d(1) = 7
          edge_d(2) = 8
          edge_d(3) = 4
       case (3)
          corner_d(0) = 3
          corner_d(1) = 2
          corner_d(2) = 6
          corner_d(3) = 7
          edge_d(0) = 2
          edge_d(1) = 7
          edge_d(2) = 9
          edge_d(3) = 10
       case (4)
          corner_d(0) = 0
          corner_d(1) = 3
          corner_d(2) = 7
          corner_d(3) = 4
          edge_d(0) = 3
          edge_d(1) = 10
          edge_d(2) = 11
          edge_d(3) = 6
       case (5)
          corner_d(0) = 4
          corner_d(1) = 5
          corner_d(2) = 6
          corner_d(3) = 7
          edge_d(0) = 5
          edge_d(1) = 8
          edge_d(2) = 9
          edge_d(3) = 11
   end select 


end subroutine caract_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine which_edge(corner,edge)

implicit none

integer, dimension(0:3), intent(IN) :: corner
integer, dimension(0:3), intent(INOUT) :: edge

integer :: i,j1,j2


do i = 0,3
  j1 = corner(i)
  if (i==3) then
    j2 = corner(0)
  else
    j2 = corner(i+1)
  endif 
  if ( (j1==0 .and. j2==1) .or. (j1==1 .and. j2==0) ) then
    edge(i) = 0  
  else if ( (j1==1 .and. j2==2) .or. (j1==2 .and. j2==1) ) then
    edge(i) = 1 
  else if ( (j1==2 .and. j2==3) .or. (j1==3 .and. j2==2) ) then
    edge(i) = 2 
  else if ( (j1==0 .and. j2==3) .or. (j1==3 .and. j2==0) ) then
    edge(i) = 3 
  else if ( (j1==1 .and. j2==5) .or. (j1==5 .and. j2==1) ) then
    edge(i) = 4 
  else if ( (j1==4 .and. j2==5) .or. (j1==5 .and. j2==4) ) then
    edge(i) = 5 
  else if ( (j1==0 .and. j2==4) .or. (j1==4 .and. j2==0) ) then
    edge(i) = 6 
  else if ( (j1==2 .and. j2==6) .or. (j1==6 .and. j2==2) ) then
    edge(i) = 7 
  else if ( (j1==5 .and. j2==6) .or. (j1==6 .and. j2==5) ) then
    edge(i) = 8 
  else if ( (j1==6 .and. j2==7) .or. (j1==7 .and. j2==6) ) then
    edge(i) = 9 
  else if ( (j1==3 .and. j2==7) .or. (j1==7 .and. j2==3) ) then
    edge(i) = 10
  else if ( (j1==4 .and. j2==7) .or. (j1==7 .and. j2==4) ) then
    edge(i) = 11
  endif 
enddo

end subroutine which_edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

