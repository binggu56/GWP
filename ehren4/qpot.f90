      subroutine qpot(nb,dv,ddv,c,am,p,q,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: nb
      real*8,intent(IN) :: w(nb),q(nb),am,p(nb),dv(nb),ddv(nb)
      complex*16,intent(in) :: c(nb)
      real*8, intent(OUT) :: u(nb),du(nb)
      complex*16 s,z0,z,psi(np)
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx,np

      u = 0d0

!      do i=1,nb
!        s = (0d0,0d0)
!        do j=1,nb
!          do k=1,nb
!            dq = q(j) - q(k)
!            z = dq/2d0-im*(p(j)+p(k))/2d0/al
!            z0 = exp(-al*z*conjg(z)+im/2d0*(p(j)+p(k))*dq)          
!            s = s+conjg(c(j))*c(k)*(z0*z+(q(k)-q(i))*z0)
!          enddo
!        enddo
!        print *,s
!        du(i) = real(s)*ddv(i)
!      enddo


!--------numerical integration---------------
      a = 1.02d0
      x0 = 1.4d0

!      do i=1,np
!        xi = xmin+dx*(i-1)
!        psi(i) = (0d0,0d0)
!        do j=1,nb
!          psi(i) = psi(i)+c(j)*exp(-al/2d0*(xi-q(j))**2+
!     +             im*p(j)*(xi-q(j)))*
!     +            dsqrt(dsqrt(al/pi))
!        enddo
!      enddo

!      do k=1,nb
!        d = 0d0
!        do i=1,np
!          xi = xmin+dx*(i-1)
!          d = d+abs(psi(i))**2*dx*2d0*(1d0-exp(-a*(xi-x0)))*a*
!     +        exp(-a*(xi-x0))
!        enddo
!        du(k) = -d
!      enddo

      d = 0d0
      do i=1,nb
        d = d+dv(i)
      enddo
      d = d/nb

      epi = 0.9d0

      n = 0d0
      x1 = 0d0
      x2 = 0d0
      do i=1,nb
        an=an+w(i)
        x1=x1+w(i)*q(i)
        x2=x2+w(i)*q(i)**2
      enddo

      b1=-0.5d0*an*an/(an*x2-x1*x1)
      b2 = 0.5d0*x1*an/(an*x2-x1*x1)
      do i=1,Nb
        r=b1*q(i)+b2
        u(i) = -(r*r+b1)/2d0/am
        du(i) = -b1*r/am-epi*dv(i)
      enddo

      return
      end subroutine
