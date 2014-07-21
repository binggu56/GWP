      subroutine qpot(Nb,ddv,c,am,p,x,w,u,du) 
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Nb
      real*8,intent(IN) :: w(Nb),x(Nb),am,p(nb),ddv(nb)
      complex*16 c(nb)
      real*8, intent(OUT) :: du(nb),u(Nb)

      an = 0d0
      x1 = 0d0
      x2 = 0d0
      do i=1,nb
        an=an+w(i)
        x1=x1+w(i)*x(i)
        x2=x2+w(i)*x(i)**2
      enddo

      b1=-0.5d0*an*an/(an*x2-x1*x1)
      b2 = 0.5d0*x1*an/(an*x2-x1*x1)
      do i=1,Nb
        r=b1*x(i)+b2
        u(i) = -(r*r+b1)/2d0/am
        du(i) = -b1*r/am
      enddo
      
      return
      end subroutine
