subroutine savetrace (Tdomain,rg,k)

! Modified by Paul Cupillard 20/09/2005


use sdomain

implicit none

type (domain), intent (IN):: Tdomain
integer, intent (IN) :: rg, k

integer :: n, ntime, nt1, nt2, ndt2
real :: dt1, dt2
character (len=100) :: fnamef


nt1 = Tdomain%TimeD%ntrace
dt1 = Tdomain%TimeD%dtmin
ndt2 = 75
nt2 = Tdomain%TimeD%ntrace/ndt2
dt2 = Tdomain%TimeD%dtmin*ndt2



do n = 0, Tdomain%n_receivers-1
      if (rg == Tdomain%sReceiver(n)%proc) then

            write (fnamef,"(a,I4.4,a,I3.3)") "STraceX/trace",n,"X",k

         open (61,file=fnamef,status="unknown",form="formatted")

            write (fnamef,"(a,I4.4,a,I3.3)") "STraceY/trace",n,"Y",k

         open (62,file=fnamef,status="unknown",form="formatted")

            write (fnamef,"(a,I4.4,a,I3.3)") "STraceZ/trace",n,"Z",k

         open (63,file=fnamef,status="unknown",form="formatted")
         if ( Tdomain%sReceiver(n)%flag == 1 ) then
          do ntime = 0, nt1-1
            write (61,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,0)
            write (62,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,1)
            write (63,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,2)
          enddo
         else if ( Tdomain%sReceiver(n)%flag==2 ) then
          do ntime = 0, nt2-1
            write (61,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,0)
            write (62,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,1)
            write (63,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,2)
          enddo
         endif
      endif
close(61); close (62); close (63)
enddo




return
end subroutine
