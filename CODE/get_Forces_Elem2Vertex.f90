subroutine get_Forces_Elem2Vertex (Tdomain,n,rg)

! modified by Elise Delavaud 26/01/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n,rg

integer :: ngllx,nglly,ngllz,i,nv
 

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
!if (nv==7) print*,'vertex 2',n,i
        select case (i)
         case (0)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,0,0,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-0', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,0,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,0,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,0,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,0,0,0:2)
            endif
         case (1)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,0,0,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-1', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,0,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,0,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,0,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,0,0,0:2)
            endif
         case (2)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-2', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,0,0:2)
            endif
         case (3)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,nglly-1,0,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-3', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,0,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,nglly-1,0,0:2)
            endif
         case (4)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,0,ngllz-1,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-4', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,0,ngllz-1,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,0,ngllz-1,0:2)
            endif
         case (5)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-5', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,0,ngllz-1,0:2)
            endif
         case (6)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-6', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,ngllz-1,0:2)
            endif
         case (7)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,0:2)
!if (rg==11 .and. nv==419) write(30,*) 'SO-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
!if (rg==11 .and. nv==3175) write(30,*) 'Neu-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
!if (rg==13 .and. nv==584) write(31,*) 'SO-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
!if (rg==13 .and. nv==3131) write(31,*) 'Neu-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
!if (rg==9 .and. nv==346) write(30,*) 'SO-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
!if (rg==9 .and. nv==3272) write(30,*) 'Neu-7', n,Tdomain%sVertex(nv)%Forces(2),Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,nglly-1,ngllz-1,0:2)
            endif
        end select
    enddo


return
end subroutine get_Forces_Elem2Vertex
