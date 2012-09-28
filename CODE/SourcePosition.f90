subroutine SourcePosition(Tdomain,rg,nb_procs)

  ! Written by Paul Cupillard 06/06/2005


  use sdomain
  use pig

  implicit none

  include 'mpif.h'

  type (domain), intent(inout) :: Tdomain
  integer, intent(in) :: rg,nb_procs

  integer :: n_src, i,j,k, n_around_elem, ix,iy,iz, ngllx,nglly,ngllz, code,n_sharing_source,nearestGLL_inno,ind_nearestGLL_inno,mat
  integer , dimension (0:20) :: el_around_node
  integer , dimension (:), allocatable :: near_node
  real :: xs,ys,zs, d,dmin, R,Xtemp,Ytemp,min_dmin,x0,y0,z0,x1,y3,z4,xi,eta,zeta,wxi,weta,wzeta,dwdxi,dwdeta,dwdzeta
!!$  double precision, dimension(0:nb_procs-1) :: distance
  real,dimension(:),allocatable::xi_vector,x_vector
  logical::into




  allocate (near_node(0:Tdomain%n_source-1))
  allocate(xi_vector(0:2))
  allocate(x_vector(0:2))
  do n_src = 0, Tdomain%n_source-1
     !Get the source coordonnates 
     if (Tdomain%curve) then
        R = Tdomain%Ssource(n_src)%Zsource
        Xtemp = tan(pi*Tdomain%Ssource(n_src)%Xsource/180)
        Ytemp = tan(pi*Tdomain%Ssource(n_src)%Ysource/180)
        D = sqrt(1 + Ytemp**2 + Xtemp**2)
        xs = R*Xtemp/D
        ys = R*Ytemp/D
        zs = R/D
     else
        xs = Tdomain%Ssource(n_src)%Xsource
        ys = Tdomain%Ssource(n_src)%Ysource
        zs = Tdomain%Ssource(n_src)%Zsource
     endif
     !Search for the elements containing the nearest nodes  to the source 
     dmin = 100000
     do i = 0,Tdomain%n_glob_nodes-1
        d = sqrt ((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
        if (d <= dmin) then
           dmin = d
           near_node(n_src) = i
        endif
     enddo
     n_around_elem = 0
     do i = 0,Tdomain%n_elem-1
        do j = 0,Tdomain%n_nodes-1
           if (Tdomain%specel(i)%Control_nodes(j)==near_node(n_src)) then
              el_around_node(n_around_elem) = i
              n_around_elem = n_around_elem + 1
           endif
        enddo
     enddo

!!$print *,'SourcePose:check1'
     !Search the minimum distance (source - GLL points)
     do i = 0,n_around_elem-1
        ngllx = Tdomain%specel(el_around_node(i))%ngllx
        nglly = Tdomain%specel(el_around_node(i))%nglly
        ngllz = Tdomain%specel(el_around_node(i))%ngllz
!!$     !allocate(Tdomain%specel(el_around_node(i))%signus(0:2))
        do ix = 0,ngllx-1
           do iy = 0,nglly-1
              do iz = 0,ngllz-1
                 j = Tdomain%specel(el_around_node(i))%Iglobnum(ix,iy,iz)
                 d = sqrt ((Tdomain%Globcoord(0,j)-xs)**2 + (Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
                 if (d < dmin) then
                    dmin = d
                 endif
              enddo
           enddo
        enddo
     enddo

!!$print *, 'SourcePos:check2'
     !Initialize the number of element sharing the source
     nearestGLL_inno=0
     !Communication to find the proc treating the subdomain containning the source

     call MPI_ALLREDUCE(dmin,min_dmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,code)
!!$print *, 'SourcePos::check3',dmin, min_dmin,rg

     if (dmin==min_dmin) then
!!$print *, 'SourcePos::check4',rg
        Tdomain%sSource(n_src)%proc(rg)%source_affected=.true.
        allocate(Tdomain%sSource(n_src)%proc(rg)%nr) 
        allocate(Tdomain%sSource(n_src)%proc(rg)%Elem(0:7))!There are at maximum 8 elements sharing a GLL-point
        print ('(" SOURCEPOSITION::Me, the ",i4,"th source, will be dealed by the proc ",i2)'),n_src+1,rg

!!$     call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,mpi_comm_world,code)
!!$     do i = 0,nb_procs-1
!!$        if (distance(i) <= dmin) then
!!$           dmin = distance(i)
!!$           Tdomain%sSource(n_src)%proc = i
!!$        endif
!!$     enddo
        !Location of the source in the element local cordinates systems(xi,eta,zeta)
        do i = 0,n_around_elem-1
           ngllx = Tdomain%specel(el_around_node(i))%ngllx
           nglly = Tdomain%specel(el_around_node(i))%nglly
           ngllz = Tdomain%specel(el_around_node(i))%ngllz
!!$     !allocate(Tdomain%specel(el_around_node(i))%signus(0:2))
           do ix = 0,ngllx-1
              do iy = 0,nglly-1
                 do iz = 0,ngllz-1
                    j = Tdomain%specel(el_around_node(i))%Iglobnum(ix,iy,iz)
                    d = sqrt ((Tdomain%Globcoord(0,j)-xs)**2 + (Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
                    if (d == dmin) then !The nearest GLL-points
                       x0=Tdomain%Coord_nodes(0,Tdomain%specel(el_around_node(i))%Control_Nodes(0))
                       x1=Tdomain%Coord_nodes(0,Tdomain%specel(el_around_node(i))%Control_Nodes(1))

                       y0=Tdomain%Coord_nodes(1,Tdomain%specel(el_around_node(i))%Control_Nodes(0))
                       y3=Tdomain%Coord_nodes(1,Tdomain%specel(el_around_node(i))%Control_Nodes(3))

                       z0=Tdomain%Coord_nodes(2,Tdomain%specel(el_around_node(i))%Control_Nodes(0))
                       z4=Tdomain%Coord_nodes(2,Tdomain%specel(el_around_node(i))%Control_Nodes(4))
!!$
!!$                       into=(xs .ge. x0) .and. (xs .le. x1) .and. &  !This comparaison is respect to the hexaedron - numerotation law of METIS 
!!$                            (ys .ge. y0) .and. (ys .le. y3) .and. &  !Verify when applicate to none-structured mesh 
!!$                            (zs .ge. z0) .and. (zs .le. z4)

                       into=(((xs-x0)*(x1-xs) .gt. -1e-12) .and.  &  !This comparaison is respect to the hexaedron - numerotation law of METIS 
                            ((ys-y0)*(y3-ys) .gt. -1e-12) .and. &  !Verify when applicate to none-structured mesh 
                            ((zs-z0)*(z4-zs) .gt. -1e-12))

!!$into= ((xs-x0)*(x1-xs) .gt. -1e-12) .and. ((ys-y0)*(y3-ys) .gt. -1e-12) .and. ((zs-z0)*(z4-zs) .gt. -1e-12)
!!$                    


!!$if (rg==4) print *, x0,x1,y0,y3,z0,z4
!if (((zs-z0)*(z4-xs) .gt. -1e-14)) 
!!$print *,'SourcePosi', ((xs-x0)*(x1-xs) .gt. -1e-12),((ys-y0)*(y3-ys) .gt. -1e-12),(((zs-z0)*(z4-zs)) .gt. -1e-12), zs,z0,z4, into
                       if (into) then ! Elements containing the source
!!$print *,'SourcePosi', x0,x1,y0,y3,(-1000-z0)*(z4+1000),xs,ys,zs,rg
                          allocate(Tdomain%sSource(n_src)%proc(rg)%Elem(nearestGLL_inno)%indspecel)
                          Tdomain%sSource(n_src)%proc(rg)%Elem(nearestGLL_inno)%indspecel=el_around_node(i)
!!$                          Tdomain%Ssource(n_src)%elem(nearestGLL_inno) = el_around_node(i)
!!$                       Tdomain%Ssource(n_src)%elem(nearestGLL_inno)%indGLLx=ix
!!$                       Tdomain%Ssource(n_src)%elem(nearestGLL_inno)%indGLLy=iy
!!$                       Tdomain%Ssource(n_src)%elem(nearestGLL_inno)%indGLLz=iz  
                          mat=Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(nearestGLL_inno)%indspecel)%mat_index 
                          xi_vector(0)=Tdomain%sSubdomain(mat)%GLLcx(ix)
                          xi_vector(1)=Tdomain%sSubdomain(mat)%GLLcx(iy)
                          xi_vector(2)=Tdomain%sSubdomain(mat)%GLLcx(iz)
                          x_vector(0)=xs-Tdomain%Globcoord(0,j)
                          x_vector(1)=ys-Tdomain%Globcoord(1,j)
                          x_vector(2)=zs-Tdomain%Globcoord(2,j)
                          call DGEMV ('N',3,3,1.,Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(nearestGLL_inno)%indspecel)%&
                               InvGrad (ix,iy,iz, 0:2,0:2),3,x_vector,1,1.,xi_vector,1)
                          Tdomain%sSource(n_src)%proc(rg)%Elem(nearestGLL_inno)%xi = xi_vector(0)
                          Tdomain%sSource(n_src)%proc(rg)%Elem(nearestGLL_inno)%eta = xi_vector(1)
                          Tdomain%sSource(n_src)%proc(rg)%Elem(nearestGLL_inno)%zeta = xi_vector(2)
                          nearestGLL_inno=nearestGLL_inno+1
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
        Tdomain%sSource(n_src)%proc(rg)%nr=nearestGLL_inno
        print ('(" SOURCEPOSITION::Always Me, the",i2, "th source, I''m now already in the proc ",i2, " and being shared by ",i2," elements.")'), n_src+1,rg,nearestGLL_inno
     endif
     call MPI_ALLREDUCE(nearestGLL_inno,n_sharing_source,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,code)    
     print *,'SOURCEPOSITION::elts sharing total',n_sharing_source,rg
     if (Tdomain%sSource(n_src)%proc(rg)%source_affected) then 
!!$        do ind_nearestGLL_inno=0,nearestGLL_inno-1
        do ind_nearestGLL_inno=0,Tdomain%sSource(n_src)%proc(rg)%nr-1
           ngllx = Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%ngllx
           nglly = Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%nglly
           ngllz = Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%ngllz
           mat=Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%mat_index
           if (Tdomain%sSource(n_src)%i_type_source == 1) then  ! Pulse directional force
              allocate  (Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%ExtForce(0:ngllx-1,0:nglly-1,0:ngllz-1))

              do iz = 0,ngllz-1
                 call pol_lagrange ('v',ngllz, Tdomain%sSubdomain(mat)%GLLcz, iz, Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%zeta,wzeta)
                 do iy = 0,nglly-1
                    call pol_lagrange ('v',nglly, Tdomain%sSubdomain(mat)%GLLcy, iy, Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%eta, weta )
                    do ix = 0,ngllx-1
                       call pol_lagrange ('v',ngllx, Tdomain%sSubdomain(mat)%GLLcx, ix, Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%xi, wxi )
                       Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%ExtForce (ix,iy,iz) = wxi*weta*wzeta/n_sharing_source
                    enddo
                 enddo
              enddo


!!$              do iz = 0,ngllz-1
!!$                 do iy = 0,nglly-1
!!$                    do ix = 0,ngllx-1
!!$                       Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%ExtForce (ix,iy,iz) = Tdomain%sSubdomain(mat)%GLLwx(ix)*&
!!$                            Tdomain%sSubdomain(mat)%GLLwy(iy)*&
!!$                            Tdomain%sSubdomain(mat)%GLLwz(iz)/n_sharing_source
!!$                    enddo
!!$                 enddo
!!$              enddo


           elseif (Tdomain%sSource(n_src)%i_type_source == 2) then !Explosion
              allocate  (Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%Explosion(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
!!$              do iz = 0,ngllz-1
!!$                 call pol_lagrange ('v',ngllz, Tdomain%sSubdomain(mat)%GLLcz, iz, zeta,wzeta)
!!$                 call pol_lagrange ('d',ngllz, Tdomain%sSubdomain(mat)%GLLcz, iz, zeta,dwdzeta)
!!$                 do iy = 0,nglly-1
!!$                    call pol_lagrange ('v',nglly, Tdomain%sSubdomain(mat)%GLLcy, iy, eta,weta)
!!$                    call pol_lagrange ('d',nglly, Tdomain%sSubdomain(mat)%GLLcy, iy, eta,dwdeta)
!!$                    do ix = 0,ngllx-1
!!$                       call pol_lagrange ('v',ngllx, Tdomain%sSubdomain(mat)%GLLcx, ix, xi, wxi )
!!$                       call pol_lagrange ('d',ngllx, Tdomain%sSubdomain(mat)%GLLcx, ix, xi, dwdxi )
!!$                       xi_vector(0)= dwdxi * weta   * wzeta
!!$                       xi_vector(1)= wxi   * dwdeta * wzeta
!!$                       xi_vector(2)= wxi   * weta   * dwdzeta 
!!$                       !DGEMV :: Y:=Alpha*A*X+Beta*Y
!!$                       call DGEMV ('N',3,3,1./n_sharing_source,&!TransA,size(A,1),size(A,2),Alpha
!!$                            Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%InvGrad (ix,iy,iz, 0:2,0:2),3,&!A,LdA
!!$                            xi_vector,1,0.,&!X,X's ind increment,Beta 
!!$                            Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%Explosion (ix,iy,iz,:),1&! Y,Y's ind increment 
!!$                            )
!!$                    enddo
!!$                 enddo
!!$              enddo

              do iz = 0,ngllz-1
                 do iy = 0,nglly-1
                    do ix = 0,ngllx-1
                       xi_vector(0)= Tdomain%sSubdomain(mat)%hprimex(ix,ix) * Tdomain%sSubdomain(mat)%GLLwy(iy) * Tdomain%sSubdomain(mat)%GLLwz(iz)
                       xi_vector(1)= Tdomain%sSubdomain(mat)%hprimey(iy,iy) * Tdomain%sSubdomain(mat)%GLLwx(ix) * Tdomain%sSubdomain(mat)%GLLwz(iz)
                       xi_vector(2)= Tdomain%sSubdomain(mat)%hprimez(iz,iz) * Tdomain%sSubdomain(mat)%GLLwy(iy) * Tdomain%sSubdomain(mat)%GLLwx(ix)
                       !DGEMV :: Y:=Alpha*A*X+Beta*Y
                       call DGEMV ('N',3,3,1./n_sharing_source,&!TransA,size(A,1),size(A,2),Alpha
                            Tdomain%specel(Tdomain%Ssource(n_src)%proc(rg)%elem(ind_nearestGLL_inno)%indspecel)%InvGrad (ix,iy,iz, 0:2,0:2),3,&!A,LdA
                            xi_vector,1,0.,&!X,X's ind increment,Beta 
                            Tdomain%sSource(n_src)%proc(rg)%Elem(ind_nearestGLL_inno)%Explosion (ix,iy,iz,:),1&! Y,Y's ind increment 
                            )
                    enddo
                 enddo
              enddo

           endif
        enddo
     endif


  enddo

  deallocate (near_node,xi_vector,x_vector)

  return
end subroutine SourcePosition

!=====================================================================
! Subroutine using the Legendre interpolation to compute the value or 
! the first derival distribued to a 1D GLL points
!=====================================================================
subroutine pol_lagrange (opt,p,GLL_grid,GLL_pos,val_in,val_out)
  implicit none

  !In-out variable 
  character(1)::opt
  integer, intent (IN) :: p,GLL_pos
  real, dimension  (0:p-1), intent (IN)::GLL_grid
  real,intent (IN)::val_in
  real, intent (OUT)::val_out

  !Local variable  
  integer :: ind0,ind1
  real :: y0


  if (p ==0 ) then
     write (*,*) "Bad number n = 0. It should be a constant!"
     stop
  endif

  select case (opt)
  case('V','v') !Value
     val_out = 1.
     if (p ==1 ) return
     do ind0 =0,p-1
        if (ind0 /= GLL_pos)  val_out = val_out * (val_in-GLL_grid(ind0))/(GLL_grid(GLL_pos)-GLL_grid(ind0))
     enddo
     return
  case ('D','d') !Derival
     val_out = 0
     do ind1 = 0 ,p-1
        y0 = 1.
        if (ind1 /= GLL_pos) then
           do ind0 = 0,p-1
              if (ind0 /= ind1 .and. ind0 /= GLL_pos)   y0 = y0 * (val_in-GLL_grid(ind0))/(GLL_grid(GLL_pos)-GLL_grid(ind0))
           enddo
           y0 = y0 / (GLL_grid(GLL_pos)-GLL_grid(ind1) ) 
           val_out = val_out + y0
        endif
     enddo
     return
  end select
end subroutine pol_lagrange

