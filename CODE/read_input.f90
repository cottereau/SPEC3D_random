subroutine read_input (Tdomain, rg)

  ! Modified by Gaetano Festa 18/02/2005
  ! Modified by Paul Cupillard 09/12/2005
  ! Modified by Elise Delavaud 15/02/2006
  ! Modified by Regis 04/2007


  use sdomain
  !use MPI
  implicit none

  include 'mpif.h'


  type (domain), intent (INOUT) :: Tdomain
  integer, intent(IN) :: rg

  logical :: logic_scheme
  logical, dimension(:), allocatable :: L_Face, L_Edge
 
  integer :: length, i, j,npml, n_aus, mat, ok, nf, ne, nv, k, nrand,nani,NRatio,ns
  real :: dtmin
  character*20 :: fnamef
  logical :: existe
  logical, dimension(:),allocatable :: FOrder


  open (11,file="input.spec",form="formatted",status="old")
  !print *,'open input.spec',rg 
  read (11,*) Tdomain%Title_simulation !== 1
  !print *,'read title',rg, Tdomain%Title_simulation
  read (11,*) Tdomain%TimeD%acceleration_scheme !== 2
  read (11,*) Tdomain%TimeD%velocity_scheme!== 3 
  read (11,*) Tdomain%TimeD%duration !== 4
  read (11,*) Tdomain%TimeD%alpha !== 5
  read (11,*) Tdomain%TimeD%beta !== 6
  read (11,*) Tdomain%TimeD%gamma !== 7
  read (11,*) Tdomain%mesh_file !== 8
  length = len_trim(Tdomain%mesh_file) + 1
  write (Tdomain%mesh_file(length:length+3),'(i4.4)') rg 
  read (11,*) Tdomain%material_file !== 9
  read (11,*) Tdomain%logicD%save_trace !== 10
  read (11,*) Tdomain%logicD%save_snapshots !== 11
  read (11,*) Tdomain%logicD%save_energy !== 12
  read (11,*) Tdomain%logicD%save_restart !== 13
  read (11,*) Tdomain%logicD%plot_grid !== 14
  read (11,*) Tdomain%logicD%run_exec !== 15
  read (11,*) Tdomain%logicD%run_debug !== 16
  read (11,*) Tdomain%logicD%run_echo !== 17
  read (11,*) Tdomain%logicD%run_restart !== 18

  if (Tdomain%logicD%save_restart) then
     read (11,*) Tdomain%TimeD%ncheck
  else 
     read( 11,*)
  endif

  if (Tdomain%logicD%save_trace) then
     read (11,*) Tdomain%station_file
     read (11,*) Tdomain%TimeD%ntrace
  else 
     read (11,*)
     read (11,*)
  endif

  if (Tdomain%logicD%save_snapshots) then
     allocate (Tdomain%Field_Order(1:7))

     Tdomain%Field_Order=.false.   
     read (11,*) Tdomain%TimeD%time_snapshots_tic,Tdomain%TimeD%time_snapshots_toc,Tdomain%TimeD%time_snapshots_step
     !print *, Tdomain%TimeD%time_snapshots_tic,Tdomain%TimeD%time_snapshots_toc,Tdomain%TimeD%time_snapshots_step,rg
     ! Introduction of fields order-form
     inquire (file="field_order_annex.in",EXIST=existe)
     if ( existe ) then  
        allocate (FOrder(1))
        open (15,file="field_order_annex.in", status="old")
        read (15,*) FOrder 
        Tdomain%Field_Order(1:1)=FOrder(1:1) !Veloc
        !print *, 'Veloc', FOrder, Tdomain%Field_Order(1)
        read (15,*) FOrder 
        Tdomain%Field_Order(2:2)=FOrder(1:1) !Stress
        !print *, 'Stress', FOrder, Tdomain%Field_Order(2)
        read (15,*) FOrder 
        Tdomain%Field_Order(3:3)=FOrder(1:1) !Energy
        !print *, 'Energy', FOrder, Tdomain%Field_Order(3)
        read (15,*), Tdomain%PML_GLL_out
        read (15,*) FOrder 
        Tdomain%Field_Order(4:4)=FOrder(1:1) !Traveltime
        read (15,*) NRatio
        allocate (Tdomain%Ratios(1:NRatio))
        read (15,*) Tdomain%Ratios
        !if (rg==0) print *,Tdomain%Ratios, Tdomain%Field_Order
        read(15,*) FOrder
        Tdomain%Field_Order(5:5)=FOrder(1:1)
        read(15,*) FOrder
        Tdomain%Field_Order(6:6)=FOrder(1:1)
        read(15,*) FOrder
        Tdomain%Field_Order(7:7)=FOrder(1:1)

        read (15,*) Tdomain%nbSlices

        if (Tdomain%Field_Order(6) .or. Tdomain%Field_Order(7) ) then
           allocate(Tdomain%Slices(Tdomain%nbSlices,6))
           Tdomain%MaxiSliceSize=5000
           allocate(Tdomain%I_elem_inSlices(Tdomain%nbSlices,Tdomain%MaxiSliceSize)) !Slice size is limited to 5000000 squares
           Tdomain%I_elem_inSlices=0.
           do ns=1,Tdomain%nbSlices
              read(15,*) Tdomain%Slices(ns,:)
           enddo
	   
        else
           do ns=1,Tdomain%nbSlices
              read(15,*) 
           enddo
        endif

        close (15)
        deallocate(FOrder)
     else
        write (*,*) "File field_order_annex.in not found, order by default only the field of velocity"
        Tdomain%Field_Order(1:1)=.true.
     endif
  else 
     read (11,*)
  endif

  if (Tdomain%logicD%save_energy) then
     read (11,*) Tdomain%TimeD%time_energy
     read (11,*) Tdomain%EnerMat
  else 
     read (11,*)
     read (11,*)
  endif

  logic_scheme = Tdomain%TimeD%acceleration_scheme .xor. Tdomain%TimeD%velocity_scheme
  if (.not. logic_scheme) then
     write (*,*) "No compatible acceleration and velocity schemes" 
     stop
  endif

  read (11,*) Tdomain%logicD%super_object
  if (Tdomain%logicD%super_object) then
     read (11,*) Tdomain%Super_object_type, Tdomain%super_object_file
  else
     read (11,*)
  endif
  read (11,*) Tdomain%logicD%Neumann
  if ( Tdomain%logicD%Neumann ) then
     read (11,*) Tdomain%neumann_file
  else
     read (11,*)
  endif


  read(11,*) Tdomain%logicD%any_source
  if (Tdomain%logicD%any_source) then
     read (11,*) Tdomain%n_source 
     !print *,'read_input::Nsource',Tdomain%n_source
     allocate (Tdomain%Ssource(0:Tdomain%n_source-1))
     if (Tdomain%n_source .gt. 5) then
!!$print *,'dqsdqsdqsdqsdqsdqsdqsdqdqsdd'
        read(11,*) !input.spec : XYZ of source
        read(11,*) !input.spec : type_source
        read(11,*) !input.spec : i_dir
        read (11,*) Tdomain%Ssource(0)%i_time_function
        read (11,*) Tdomain%Ssource(0)%tau_b
        read (11,*) Tdomain%Ssource(0)%cutoff_freq
        open (111,file="sources.posi", form="formatted", status="unknown")
        read (111,*) ! Legends
        do i = 0, Tdomain%n_source - 1
           read (111,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource,Tdomain%Ssource(i)%i_type_source,Tdomain%Ssource(i)%i_dir   
           Tdomain%Ssource(i)%i_time_function=Tdomain%Ssource(0)%i_time_function
           Tdomain%Ssource(i)%tau_b= Tdomain%Ssource(0)%tau_b
           Tdomain%Ssource(i)%cutoff_freq= Tdomain%Ssource(0)%cutoff_freq
        enddo
        close(111)
     else
        do i = 0, Tdomain%n_source - 1
           read (11,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource
           read (11,*) Tdomain%Ssource(i)%i_type_source
           if (Tdomain%Ssource(i)%i_type_source == 1 ) then
              read (11,*) Tdomain%Ssource(i)%i_dir
           else
              read (11,*) 
              Tdomain%Ssource(i)%i_dir=4 
           endif
!!$print *, Tdomain%Ssource(i)%i_type_source,'i_type_source'
           read (11,*) Tdomain%Ssource(i)%i_time_function
           read (11,*) Tdomain%Ssource(i)%tau_b
           read (11,*) Tdomain%Ssource(i)%cutoff_freq
        enddo
     endif
  else 
     read(11,*)  !input.spec : number of sources
     read(11,*)  !input.spec : XYZ of source
     read(11,*)  !input.spec : type_source
     read(11,*)  !input.spec : i_dir
     read (11,*) !input.spec : time_function
     read (11,*) !input.spec : tau_b
     read (11,*) !input.spec : cutoff_freq
  endif
     read (11,*) Tdomain%SaveMecaField 
     print *, 'Read_input:SaveMecaField',Tdomain%SaveMecaField 
     read (11,*) Tdomain%TimeReverse
     if (Tdomain%TimeReverse) read (11,*) Tdomain%TimeD%Ntime2reverse
  close (11)


  ! If echo modality write the read parameter in a file
  if (Tdomain%logicD%run_echo) then
     open (91,file="input_spec_echo", form="formatted", status="unknown")
     write (91,*) Tdomain%Title_simulation
     write (91,*) Tdomain%TimeD%acceleration_scheme
     write (91,*) Tdomain%TimeD%velocity_scheme
     write (91,*) Tdomain%TimeD%duration
     write (91,*) Tdomain%TimeD%alpha
     write (91,*) Tdomain%TimeD%beta
     write (91,*) Tdomain%TimeD%gamma
     write (91,*) Tdomain%mesh_file
     write (91,*) Tdomain%material_file
     write (91,*) Tdomain%logicD%save_trace
     write (91,*) Tdomain%logicD%save_snapshots
     write (91,*) Tdomain%logicD%save_energy
     write (91,*) Tdomain%logicD%plot_grid
     write (91,*) Tdomain%logicD%run_exec
     write (91,*) Tdomain%logicD%run_debug
     write (91,*) Tdomain%logicD%run_echo

     if (Tdomain%logicD%save_trace) then
        write (91,*) Tdomain%station_file
     else 
        write (91,*) " No parameter ned here"
     endif

     if (Tdomain%logicD%save_snapshots) then
        write (91,*) Tdomain%TimeD%time_snapshots_tic,Tdomain%TimeD%time_snapshots_toc,Tdomain%TimeD%time_snapshots_step
     else 
        write (91,*) " No parameter ned here"
     endif

     



if (Tdomain%logicD%super_object) then
        write (91,*) Tdomain%Super_object_type,"   ", Tdomain%super_object_file
     else
        write (91,*) "no super objects present"
     endif

     if (Tdomain%logicD%any_source) then
        write (91,*) Tdomain%n_source 
        do i = 0, Tdomain%n_source - 1
           write (91,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource
           write (91,*) Tdomain%Ssource(i)%i_type_source
           if (Tdomain%Ssource(i)%i_type_source == 1 ) then
              write (91,*) Tdomain%Ssource(i)%i_dir
           else
              write (91,*) "No parameter need here"
           endif
           write (91,*) Tdomain%Ssource(i)%i_time_function
           write (91,*) Tdomain%Ssource(i)%tau_b
           write (91,*) Tdomain%Ssource(i)%cutoff_freq
        enddo
     else
        write (91,*)  "No available sources "
     endif
     write (91,*) "All right, runner ?"
     close (91)
  endif


  ! Read Mesh properties
  open (12, file=Tdomain%mesh_file, iostat=ok, status="old", form="formatted")
  if (ok/=0) then
     write (*,*) "Process ",rg, " can't open his mesh_file"
     stop
  endif
  read (12,*) Tdomain%n_dime
  read (12,*) Tdomain%n_glob_nodes
  read (12,*) Tdomain%curve
  allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
  do i = 0,Tdomain%n_glob_nodes-1
     read (12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
  enddo
  read (12,*) Tdomain%n_elem
  allocate (Tdomain%specel(0:Tdomain%n_elem-1))
  read (12,*) Tdomain%n_mat
  do i = 0, Tdomain%n_elem - 1
     read(12,*) Tdomain%specel(i)%mat_index
  enddo
  if (Tdomain%n_dime == 3) then
     read (12,*) Tdomain%n_nodes
     do i = 0, Tdomain%n_elem - 1
        allocate (Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        read(12,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1)
     enddo
     read (12,*) Tdomain%n_face
     allocate (Tdomain%sFace(0:Tdomain%n_face-1))
     do i = 0, Tdomain%n_elem - 1
        read(12,*) Tdomain%specel(i)%Near_Faces(0:5)
        read(12,*) Tdomain%specel(i)%Orient_Faces(0:5)
     enddo
     read (12,*) Tdomain%n_edge
     allocate (Tdomain%sEdge(0:Tdomain%n_edge-1))
     do i = 0, Tdomain%n_elem - 1
        read(12,*) Tdomain%specel(i)%Near_Edges(0:11)
        read(12,*) Tdomain%specel(i)%Orient_Edges(0:11)
     enddo
     read (12,*) Tdomain%n_vertex
     allocate (Tdomain%sVertex(0:Tdomain%n_vertex-1))
     do i = 0, Tdomain%n_elem - 1
        read(12,*) Tdomain%specel(i)%Near_Vertices(0:7)
     enddo
     read(12,*)
     ! Information about super-object
     read(12,*) ! Super Object
     read(12,*) Tdomain%logicD%super_object_local_present
     if (Tdomain%logicD%super_object) then
        if ( Tdomain%logicD%super_object_local_present ) then
           read(12,*) ! Faces in Super Object
           read(12,*) Tdomain%sPlaneW%n_faces
           allocate (Tdomain%sPlaneW%pFace(0:Tdomain%sPlaneW%n_faces-1))
           read(12,*) ! 4 Edges for each face
           do nf = 0, Tdomain%sPlaneW%n_faces-1 
              read(12,*) Tdomain%sPlaneW%pFace(nf)%Near_Edges(0:3)
              read(12,*) Tdomain%sPlaneW%pFace(nf)%Orient_Edges(0:3)
           enddo
           read(12,*) ! 4 Vertices for each face
           do nf = 0, Tdomain%sPlaneW%n_faces-1 
              read(12,*) Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0:3)
           enddo
           read(12,*) ! Glob number of up and down faces and orientation of down compared to up
           read(12,*)
           do nf = 0, Tdomain%sPlaneW%n_faces-1
              read(12,*) Tdomain%sPlaneW%pFace(nf)%Face_UP, Tdomain%sPlaneW%pFace(nf)%Face_DOWN, Tdomain%sPlaneW%pFace(nf)%Orient
           enddo
           read(12,*) ! Glob number of up and down edges and orientation of down compared to up
           read(12,*) Tdomain%sPlaneW%n_edges
           allocate (Tdomain%sPlaneW%pEdge(0:Tdomain%sPlaneW%n_edges-1))
           do ne = 0, Tdomain%sPlaneW%n_edges-1
              read(12,*) Tdomain%sPlaneW%pEdge(ne)%Edge_UP, Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN, Tdomain%sPlaneW%pEdge(ne)%Orient
           enddo
           read(12,*) ! Glob number of up and down vertices and orientation of down compared to up
           read(12,*) Tdomain%sPlaneW%n_vertices
           allocate (Tdomain%sPlaneW%pVertex(0:Tdomain%sPlaneW%n_vertices-1))
           do nv = 0, Tdomain%sPlaneW%n_vertices-1
              read(12,*) Tdomain%sPlaneW%pVertex(nv)%Vertex_UP, Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
           enddo
        endif
     endif
     read(12,*)
     read(12,*) ! Neumann
     read(12,*) Tdomain%logicD%Neumann_local_present
     if ( Tdomain%logicD%Neumann_local_present ) then
        read(12,*) ! Faces in Neumann
        read(12,*) Tdomain%sNeu%n_faces
        allocate (Tdomain%sNeu%nFace(0:Tdomain%sNeu%n_faces-1))
        read(12,*) ! 4 Edges for each face
        do nf = 0, Tdomain%sNeu%n_faces-1 
           read(12,*) Tdomain%sNeu%nFace(nf)%Near_Edges(0:3)
           read(12,*) Tdomain%sNeu%nFace(nf)%Orient_Edges(0:3)
        enddo
        read(12,*) ! 4 Vertices for each face
        do nf = 0, Tdomain%sNeu%n_faces-1 
           read(12,*) Tdomain%sNeu%nFace(nf)%Near_Vertices(0:3)
        enddo
        read(12,*) ! Glob number of faces
        read(12,*)
        do nf = 0, Tdomain%sNeu%n_faces-1
           read(12,*) Tdomain%sNeu%nFace(nf)%Face
        enddo
        read(12,*) ! Glob number of edges 
        read(12,*) Tdomain%sNeu%n_edges
        allocate (Tdomain%sNeu%nEdge(0:Tdomain%sNeu%n_edges-1))
        do ne = 0, Tdomain%sNeu%n_edges-1
           read(12,*) Tdomain%sNeu%nEdge(ne)%Edge
        enddo
        read(12,*) ! Glob number vertices 
        read(12,*) Tdomain%sNeu%n_vertices
        allocate (Tdomain%sNeu%nVertex(0:Tdomain%sNeu%n_vertices-1))
        do nv = 0, Tdomain%sNeu%n_vertices-1
           read(12,*) Tdomain%sNeu%nVertex(nv)%Vertex
        enddo
     endif
     read(12,*)
     read (12,*) Tdomain%n_proc
     allocate (Tdomain%sComm(0:Tdomain%n_proc-1))
     do i = 0,Tdomain%n_proc-1
        read(12,*) Tdomain%sComm(i)%nb_faces, Tdomain%sComm(i)%nb_edges, Tdomain%sComm(i)%nb_vertices, &
             Tdomain%sComm(i)%nb_edges_so, Tdomain%sComm(i)%nb_vertices_so, Tdomain%sComm(i)%nb_edges_neu, Tdomain%sComm(i)%nb_vertices_neu     
        if (Tdomain%sComm(i)%nb_faces>0) then
           allocate (Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
           allocate (Tdomain%sComm(i)%orient_faces(0:Tdomain%sComm(i)%nb_faces-1))
           do j = 0,Tdomain%sComm(i)%nb_faces-1
              read(12,*) Tdomain%sComm(i)%faces(j),Tdomain%sComm(i)%orient_faces(j)
           enddo
        endif
        if (Tdomain%sComm(i)%nb_edges>0) then
           allocate (Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
           allocate (Tdomain%sComm(i)%orient_edges(0:Tdomain%sComm(i)%nb_edges-1))
           do j = 0,Tdomain%sComm(i)%nb_edges-1
              read(12,*) Tdomain%sComm(i)%edges(j),Tdomain%sComm(i)%orient_edges(j)   
           enddo
        endif
        if (Tdomain%sComm(i)%nb_vertices>0) then
           allocate (Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
           do j = 0,Tdomain%sComm(i)%nb_vertices-1
              read(12,*) Tdomain%sComm(i)%vertices(j)
           enddo
        endif
        if (Tdomain%sComm(i)%nb_edges_so>0) then
           allocate (Tdomain%sComm(i)%edges_SO(0:Tdomain%sComm(i)%nb_edges_so-1))
           allocate (Tdomain%sComm(i)%orient_edges_SO(0:Tdomain%sComm(i)%nb_edges_so-1))
           do j = 0,Tdomain%sComm(i)%nb_edges_so-1
              read(12,*) Tdomain%sComm(i)%edges_SO(j),Tdomain%sComm(i)%orient_edges_SO(j)          
           enddo
        endif
        if (Tdomain%sComm(i)%nb_vertices_so>0) then
           allocate (Tdomain%sComm(i)%vertices_SO(0:Tdomain%sComm(i)%nb_vertices_so-1))
           do j = 0,Tdomain%sComm(i)%nb_vertices_so-1
              read(12,*) Tdomain%sComm(i)%vertices_SO(j)
           enddo
        endif
        if (Tdomain%sComm(i)%nb_edges_neu>0) then
           allocate (Tdomain%sComm(i)%edges_Neu(0:Tdomain%sComm(i)%nb_edges_neu-1))
           allocate (Tdomain%sComm(i)%orient_edges_Neu(0:Tdomain%sComm(i)%nb_edges_neu-1))
           do j = 0,Tdomain%sComm(i)%nb_edges_neu-1
              read(12,*) Tdomain%sComm(i)%edges_Neu(j),Tdomain%sComm(i)%orient_edges_Neu(j)          
           enddo
        endif
        if (Tdomain%sComm(i)%nb_vertices_neu>0) then
           allocate (Tdomain%sComm(i)%vertices_Neu(0:Tdomain%sComm(i)%nb_vertices_neu-1))
           do j = 0,Tdomain%sComm(i)%nb_vertices_neu-1
              read(12,*) Tdomain%sComm(i)%vertices_Neu(j)
           enddo
        endif
     enddo
     read(12,*)
     read(12,*) ! Glob number vertices 
     read(12,*) Tdomain%logicD%Save_Surface
     if ( Tdomain%logicD%Save_Surface ) then
        read(12,*) Tdomain%sSurf%n_vertices
        allocate (Tdomain%sSurf%nVertex(0:Tdomain%sSurf%n_vertices-1))
        do nv = 0, Tdomain%sSurf%n_vertices-1
           read(12,*) Tdomain%sSurf%nVertex(nv)%Vertex
        enddo
     endif
  else
     write (*,*) "A dimension different from 3 is not yet taken into account"
     stop
  endif
  close (12)


  npml = 0
  nrand = 0
  nani=0
  allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))
  open (13, file=Tdomain%material_file, status="old", form="formatted")
  print *,'Material file ', Tdomain%material_file  
  read (13,*) n_aus
  if (n_aus /= Tdomain%n_mat) then
     write (*,*) "Incompatibility between the mesh file and the material file for n_mat "
     stop
  endif
  do i = 0,Tdomain%n_mat-1
     Tdomain%sSubDomain(i)%aniso=.false.
     read (13,*) Tdomain%sSubDomain(i)%material_type
!!$     print *,  Tdomain%sSubDomain(i)%material_type
!!$     print *,   Tdomain%sSubDomain(i)%material_type(2:2)
     if (Tdomain%sSubDomain(i)%material_type(2:2)=="I") then
        read(13,*) Tdomain%sSubDomain(i)%Pspeed, &
             Tdomain%sSubDomain(i)%Sspeed,&
             Tdomain%sSubDomain(i)%dDensity, &
             Tdomain%sSubDomain(i)%NGLLx,&
             Tdomain%sSubDomain(i)%NGLLy, &
             Tdomain%sSubDomain(i)%NGLLz, Tdomain%sSubDomain(i)%Dt
     elseif (Tdomain%sSubDomain(i)%material_type(2:2)=="A") then 
        read(13,*) Tdomain%sSubDomain(i)%Coefs21, Tdomain%sSubDomain(i)%dDensity, &
             Tdomain%sSubDomain(i)%NGLLx,&
             Tdomain%sSubDomain(i)%NGLLy, &
             Tdomain%sSubDomain(i)%NGLLz, Tdomain%sSubDomain(i)%Dt
        Tdomain%sSubDomain(i)%aniso=.true.
!!$        print *, Tdomain%sSubDomain(i)%Coefs21(1),Tdomain%sSubDomain(i)%Coefs21(3),Tdomain%sSubDomain(i)%Coefs21(15),'----',i
        nani=nani+1
     else
        write (*,*) "Material type (anisotropic or isotropic) not precised "
        stop
     endif

     if (Tdomain%sSubDomain(i)%material_type(1:1) == "P")  then
        Tdomain%sSubDomain(i)%wpml = npml
        npml = npml + 1
     elseif (Tdomain%sSubDomain(i)%material_type(1:1) == "R") then
        nrand = nrand+1
     elseif (Tdomain%sSubDomain(i)%material_type(1:1) == "M") then
        npml = npml+1
        nrand = nrand+1
     endif
  enddo

  Tdomain%any_PML = .false.
  Tdomain%any_FPML = .false.
  if (npml > 0) then
     Tdomain%any_PML = .true.
     read(13,*); read(13,*)
     do i = 0,Tdomain%n_mat-1
        if ((Tdomain%sSubdomain(i)%material_type(1:1) == "P").or. &
             (Tdomain%sSubdomain(i)%material_type(1:1) == "M")) then
           read (13,*) Tdomain%sSubdomain(i)%Filtering, Tdomain%sSubdomain(i)%npow, &
                Tdomain%sSubdomain(i)%Apow, Tdomain%sSubdomain(i)%Px, &
                Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Py, &
                Tdomain%sSubdomain(i)%Forward, Tdomain%sSubdomain(i)%Pz, &
                Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq,&
                Tdomain%sSubdomain(i)%coPMLrate
           if (Tdomain%sSubdomain(i)%Filtering) Tdomain%any_FPML = .true.
           Tdomain%sSubdomain(i)%PML_GLL_out=Tdomain%PML_GLL_out
        endif
     enddo
  endif
  Tdomain%any_random = .false.
  if (nrand > 0) then
     Tdomain%any_random = .true.
     read(13,*); read(13,*)
     do i = 0,Tdomain%n_mat-1
        Tdomain%sSubdomain(i)%random = .false.
        if ((Tdomain%sSubdomain(i)%material_type(1:1) == "R").or. &         
             (Tdomain%sSubdomain(i)%material_type(1:1) == "M")) then         
           read (13,*) Tdomain%sSubdomain(i)%random_type, &
                Tdomain%sSubdomain(i)%Correlation_type, &
                Tdomain%sSubdomain(i)%delta, &
                Tdomain%sSubdomain(i)%CorLx, &
                Tdomain%sSubdomain(i)%CorLy, &
                Tdomain%sSubdomain(i)%CorLz, &
                Tdomain%sSubdomain(i)%CorNx, &
                Tdomain%sSubdomain(i)%CorNy, &
                Tdomain%sSubdomain(i)%CorNz, &
                Tdomain%sSubdomain(i)%ChosenSeed, &               
                Tdomain%sSubdomain(i)%SpheDev, &!--------------------------------------
                Tdomain%sSubdomain(i)%delta_kernel, &
                Tdomain%sSubdomain(i)%random_type_kernel, &!--------------------------------------New_random
                Tdomain%sSubdomain(i)%reload_mat!-------------Reload_mat
           
           if(Tdomain%sSubdomain(i)%material_type(1:1) == "R") Tdomain%sSubdomain(i)%material_type(1:1) = "S"
           Tdomain%sSubdomain(i)%random = .true.
        endif
     enddo
  endif
  close (13)

  Tdomain%any_aniso = .false.
  if (nani>0) Tdomain%any_aniso = .true. 


  allocate (L_Face(0:Tdomain%n_face-1))
  L_Face = .true.
  allocate (L_Edge(0:Tdomain%n_edge-1))
  L_Edge = .true.
  do i = 0,Tdomain%n_elem-1
     mat = Tdomain%specel(i)%mat_index
     Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
     Tdomain%specel(i)%nglly = Tdomain%sSubDomain(mat)%NGLLy
     Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
     do j = 0,5
        nf = Tdomain%specel(i)%Near_Faces(j)
        if (L_Face(nf)) then
           L_Face(nf) = .false.
           if (j==0 .or. j==5) then
              Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
              Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%nglly
              Tdomain%sFace(nf)%dir = j 
           else if (j==1 .or. j==3) then
              Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
              Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
              Tdomain%sFace(nf)%dir = j
           else
              Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%nglly
              Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
              Tdomain%sFace(nf)%dir = j 
           endif
        endif
     enddo
     do j = 0,11
        ne = Tdomain%specel(i)%Near_Edges(j)
        if (L_Edge(ne)) then
           L_Edge(ne) = .false.
           if (j==0 .or. j==2 .or. j==5 .or. j==9) then
              Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllx
           else if (j==1 .or. j==3 .or. j==8 .or. j==11) then
              Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%nglly
           else
              Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllz
           endif
        endif
     enddo
     
  enddo
  deallocate (L_Face,L_Edge)

  if (Tdomain%logicD%super_object_local_present) then
     do nf = 0, Tdomain%sPlaneW%n_faces-1
        n_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
        Tdomain%sPlaneW%pFace(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
        Tdomain%sPlaneW%pFace(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
        Tdomain%sPlaneW%pFace(nf)%dir = Tdomain%sFace(n_aus)%dir
     enddo
     do ne = 0, Tdomain%sPlaneW%n_edges-1
        n_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
        Tdomain%sPlaneW%pEdge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
     enddo
  endif
  if (Tdomain%logicD%neumann_local_present) then
     do nf = 0, Tdomain%sNeu%n_faces-1
        n_aus = Tdomain%sNeu%nFace(nf)%Face
        Tdomain%sNeu%nFace(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
        Tdomain%sNeu%nFace(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
        Tdomain%sNeu%nFace(nf)%dir = Tdomain%sFace(n_aus)%dir
     enddo
     do ne = 0, Tdomain%sNeu%n_edges-1
        n_aus = Tdomain%sNeu%nEdge(ne)%Edge
        Tdomain%sNeu%nEdge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
     enddo
  endif


  dtmin = 1e20 
  do i = 0,Tdomain%n_mat-1
     if (Tdomain%sSubDomain(i)%Dt < dtmin) dtmin = Tdomain%sSubDomain(i)%Dt
  enddo
  Tdomain%TimeD%dtmin = dtmin
  if (dtmin > 0) then
     Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/dtmin)
  else
     write (*,*) "Your dt min is zero : verify it"
     stop
  endif


  do i = 0, Tdomain%n_mat-1
     if (.not. Tdomain%sSubDomain(i)%aniso) call Lame_coefficients (Tdomain%sSubDomain(i))
  enddo


  if (Tdomain%logicD%save_trace) then
     open (14, file=Tdomain%station_file, status="old")
     read (14,*) Tdomain%n_receivers
     allocate (Tdomain%sReceiver(0:Tdomain%n_receivers-1))
     do i = 0, Tdomain%n_receivers-1
        read(14,*) Tdomain%sReceiver(i)%xRec, Tdomain%sReceiver(i)%yRec, Tdomain%sReceiver(i)%zRec, Tdomain%sReceiver(i)%flag
     enddo
     close (14)
  endif

  ! PostProcess Gid
  write (fnamef,"(a,I2.2,a)") "Proc",rg,".flavia.msh"
  open (22, file=fnamef, iostat=ok, status="unknown", form="formatted")
  write(22,*) Tdomain%n_glob_nodes,Tdomain%n_elem
  do i = 0,Tdomain%n_glob_nodes-1
!!$     write (22,"(I4.4,3f)") i+1,(Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
write (22,"(I4.4,a,f,a,f,a,f)") i+1,' ',Tdomain%Coord_nodes(0,i),' ',Tdomain%Coord_nodes(1,i),' ',Tdomain%Coord_nodes(2,i)
  enddo
  write(22,*)
  do i = 0, Tdomain%n_elem - 1
     write (22,"(I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I2.2)") i+1,' ',&
          Tdomain%specel(i)%Control_Nodes(0)+1,' ',Tdomain%specel(i)%Control_Nodes(1)+1,' ',&
          Tdomain%specel(i)%Control_Nodes(2)+1,' ',Tdomain%specel(i)%Control_Nodes(3)+1,' ',&
          Tdomain%specel(i)%Control_Nodes(4)+1,' ',Tdomain%specel(i)%Control_Nodes(5)+1,' ',&
          Tdomain%specel(i)%Control_Nodes(6)+1,' ',Tdomain%specel(i)%Control_Nodes(7)+1,' ',&
          Tdomain%specel(i)%mat_index+1 
  enddo
  close(22)


  do i = 0,Tdomain%n_face-1
     do j = 0, Tdomain%n_elem - 1
        do k=0,5
           if ( Tdomain%specel(j)%Near_Faces(k) == i )  Tdomain%sFace(i)%Which_Elem = j
        enddo
     enddo
  enddo

  return
end subroutine read_Input

! $Id: read_input.f90,v 1.6 2007/08/30 15:00:58 taquanganh Exp $ !
