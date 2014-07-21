      subroutine qpot(Nb,ddv,c,am,p,x,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Nb
      real*8, intent(OUT) :: du(Nb),u(Nb)      
      integer INFO
      real*8 :: w(Nb),x(Nb),c(nb),p(nb)
      real*8 :: f(2),rx(Nb),ddv(nb),s(2,2),cr(2,1)
! define c matrix 
        cr(1,1)=-0.5d0
        cr(2,1)=0d0

! quantum force
        s=0d0
        do i=1,Ntraj
          f=(/x(i),1d0/)
        do m=1,2
          do n=1,2
            s(m,n)=w(i)*f(m)*f(n)+s(m,n)
          end do
        end do
        end do

! calculate matrix c(t)
        call DPOSV('U',2,1,s,2,cr,2,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

! the momentum operator r=cf
        do i=1,Ntraj
          rx(i)=cr(1,1)*x(i)+cr(2,1)

! calculate quantum potential
          u(i)=-rx(i)**2/(2d0*am)-cr(1,1)/(2d0*am)
          du(i)=-rx(i)*cr(1,1)/am
        end do

        return
        end subroutine
