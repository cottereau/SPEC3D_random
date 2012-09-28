subroutine deallocate_domain (Tdomain,rg)

  ! Modified by Gaetano Festa 25/02/2005
  ! Modified by Paul Cupillard 08/12/2005


  use sdomain

  implicit none

  type(domain), intent (INOUT):: Tdomain

  integer :: n,mat,rg

  nullify (Tdomain%GlobCoord)
  nullify (Tdomain%Coord_Nodes)
!print *, 'Nullify::check1',rg
  do n = 0,Tdomain%n_elem-1
     mat = Tdomain%specel(n)%mat_index
     Tdomain%specel(n)%aniso = .false.
 !print *, 'Nullify::check2',rg
     if (Tdomain%sSubDomain(mat)%material_type(2:2) == "A") Tdomain%specel(n)%aniso = .true.
 !print *, 'Nullify::check3',rg
     nullify (Tdomain%specel(n)%Density)
     nullify (Tdomain%specel(n)%MassMat)
     nullify (Tdomain%specel(n)%IglobNum)
     nullify (Tdomain%specel(n)%Control_Nodes)
     nullify (Tdomain%specel(n)%Jacob)
 !print *, 'Nullify::check4',rg
     if (Tdomain%TimeD%velocity_scheme) then
        nullify (Tdomain%specel(n)%Veloc )
        nullify (Tdomain%specel(n)%Accel)
        nullify (Tdomain%specel(n)%V0 )  
        nullify (Tdomain%specel(n)%Forces)
 !print *, 'Nullify::check5',rg
        if (Tdomain%specel(n)%PML) then
           nullify (Tdomain%specel(n)%Acoeff)
           nullify (Tdomain%specel(n)%Diagonal_Stress)
           nullify (Tdomain%specel(n)%Diagonal_Stress1)
           nullify (Tdomain%specel(n)%Diagonal_Stress2)
           nullify (Tdomain%specel(n)%Diagonal_Stress3)
           nullify (Tdomain%specel(n)%Residual_Stress)
           nullify (Tdomain%specel(n)%Residual_Stress1)
           nullify (Tdomain%specel(n)%Residual_Stress2)
 !print *, 'Nullify::check6',rg
           if  (Tdomain%specel(n)%aniso) then
              nullify (Tdomain%specel(n)%Diagonal_Stress4a)
              nullify (Tdomain%specel(n)%Diagonal_Stress4b)
              nullify (Tdomain%specel(n)%Diagonal_Stress5a)
              nullify (Tdomain%specel(n)%Diagonal_Stress5b)
              nullify (Tdomain%specel(n)%Diagonal_Stress6a)
              nullify (Tdomain%specel(n)%Diagonal_Stress6b)
              nullify (Tdomain%specel(n)%Residual_Stress3)
              nullify (Tdomain%specel(n)%Residual_Stress4a)
              nullify (Tdomain%specel(n)%Residual_Stress4b)
              nullify (Tdomain%specel(n)%Residual_Stress5a)
              nullify (Tdomain%specel(n)%Residual_Stress5b)
              nullify (Tdomain%specel(n)%Residual_Stress6a)
              nullify (Tdomain%specel(n)%Residual_Stress6b)
 !print *, 'Nullify::check7',rg
           endif
           nullify (Tdomain%specel(n)%Veloc1)
           nullify (Tdomain%specel(n)%Veloc2)
           nullify (Tdomain%specel(n)%Veloc3)
           nullify (Tdomain%specel(n)%Forces1)
           nullify (Tdomain%specel(n)%Forces2)
           nullify (Tdomain%specel(n)%Forces3)
           nullify (Tdomain%specel(n)%DumpSx)
           nullify (Tdomain%specel(n)%DumpSy)
           nullify (Tdomain%specel(n)%DumpSz)
           nullify (Tdomain%specel(n)%DumpVx)
           nullify (Tdomain%specel(n)%DumpVy)
           nullify (Tdomain%specel(n)%DumpVz)
 !print *, 'Nullify::check8',rg
        else
           nullify (Tdomain%specel(n)%Acoeff)
           nullify (Tdomain%specel(n)%Displ ) 
!!$           if (allocated(Tdomain%specel(n)%c)) nullify(Tdomain%specel(n)%c)
!!$           if (allocated(Tdomain%specel(n)%MaxAbsVeloc)) &
!!$                nullify(Tdomain%specel(n)%MaxAbsVeloc,Tdomain%specel(n)%TravelTime,Tdomain%specel(n)%TravelTimeFound) 
        endif
     endif
  enddo
! print *, 'Nullify::check9',rg
  do n = 0, Tdomain%n_face-1
     nullify (Tdomain%sFace(n)%MassMat)
     nullify (Tdomain%sFace(n)%Veloc)
     nullify (Tdomain%sFace(n)%Forces)
     nullify (Tdomain%sFace(n)%Accel)
     nullify (Tdomain%sFace(n)%V0)  
     if (Tdomain%sFace(n)%PML) then 
        nullify (Tdomain%sFace(n)%Forces1) 
        nullify (Tdomain%sFace(n)%Forces2) 
        nullify (Tdomain%sFace(n)%Forces3) 
        nullify (Tdomain%sFace(n)%Veloc1) 
        nullify (Tdomain%sFace(n)%Veloc2)
        nullify (Tdomain%sFace(n)%Veloc3)  
        nullify (Tdomain%sFace(n)%DumpVx) 
        nullify (Tdomain%sFace(n)%DumpVy) 
        nullify (Tdomain%sFace(n)%DumpVz) 
     else
        nullify (Tdomain%sFace(n)%Displ)
     endif
  enddo
 !print *, 'Nullify::check10',rg
  do n = 0,Tdomain%n_edge-1
     nullify (Tdomain%sEdge(n)%MassMat)
     nullify (Tdomain%sEdge(n)%Veloc)
     nullify (Tdomain%sEdge(n)%Forces)
     nullify (Tdomain%sEdge(n)%Accel)
     nullify (Tdomain%sEdge(n)%V0)
     if (Tdomain%sEdge(n)%PML) then
        nullify (Tdomain%sEdge(n)%Forces1)
        nullify (Tdomain%sEdge(n)%Forces2)
        nullify (Tdomain%sEdge(n)%Forces3)
        nullify (Tdomain%sEdge(n)%Veloc1)
        nullify (Tdomain%sEdge(n)%Veloc2)
        nullify (Tdomain%sEdge(n)%Veloc3)
        nullify (Tdomain%sEdge(n)%DumpVx)
        nullify (Tdomain%sEdge(n)%DumpVy)
        nullify (Tdomain%sEdge(n)%DumpVz)
     else
        nullify (Tdomain%sEdge(n)%Displ)
     endif
  enddo
 !print *, 'Nullify::check11',rg
  do n = 0,Tdomain%n_proc-1
     if (Tdomain%sComm(n)%ngll>0) then
        nullify (Tdomain%sComm(n)%GiveForces)
        nullify (Tdomain%sComm(n)%TakeForces)
     endif
     !    if (Tdomain%sComm(n)%ngllPML>0) then
     !        nullify (Tdomain%sComm(n)%GiveForcesPML)
     !        nullify (Tdomain%sComm(n)%TakeForcesPML)
     !    endif
  enddo
 !print *, 'Nullify::check12',rg
  do n = 0, Tdomain%n_mat-1 
     if (associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcx) .or. &
          associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcy )) then
        nullify (Tdomain%sSubdomain(n)%GLLcz)
        nullify(Tdomain%sSubdomain(n)%GLLwz)
        nullify(Tdomain%sSubdomain(n)%hprimez)
        nullify(Tdomain%sSubdomain(n)%hTprimez)
     else
        nullify (Tdomain%sSubdomain(n)%GLLcz)
        nullify (Tdomain%sSubdomain(n)%GLLwz)
        nullify (Tdomain%sSubdomain(n)%hprimez)
        nullify (Tdomain%sSubdomain(n)%hTprimez)
     endif
     if (associated(Tdomain%sSubdomain(n)%GLLcy, Tdomain%sSubdomain(n)%GLLcx) ) then
        nullify (Tdomain%sSubdomain(n)%GLLcy)
        nullify(Tdomain%sSubdomain(n)%GLLwy)
        nullify(Tdomain%sSubdomain(n)%hprimey)
        nullify(Tdomain%sSubdomain(n)%hTprimey)
     else
        nullify (Tdomain%sSubdomain(n)%GLLcy)
        nullify (Tdomain%sSubdomain(n)%GLLwy)
        nullify (Tdomain%sSubdomain(n)%hprimey)
        nullify (Tdomain%sSubdomain(n)%hTprimey)
     endif
     nullify (Tdomain%sSubdomain(n)%GLLcx)
     nullify (Tdomain%sSubdomain(n)%GLLwx)
     nullify (Tdomain%sSubdomain(n)%hprimex)
     nullify (Tdomain%sSubdomain(n)%hTprimex)
  enddo
!print *, 'Nullify::check13',rg
!!$  do n = 0, Tdomain%n_receivers-1
!!$     nullify (Tdomain%sReceiver(n)%StoreTrace)
!!$  enddo
!print *, 'Nullify::check14',rg
  return
end subroutine deallocate_domain
