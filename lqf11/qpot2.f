      subroutine qpot(nb,ddv,c,am,p,q,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: nb
      real*8,intent(IN) :: w(nb),q(nb),am,p(nb),ddv(nb)
      complex*16,intent(in) :: c(nb)
      real*8, intent(OUT) :: u(nb),du(nb)
      complex*16 s,z0,z,psi(np)
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx,np

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
      anrm = dsqrt(dsqrt(al/pi))
      h = 1d-3
      do i=1,nb
        psi1 = (0d0,0d0)
        psi0 = (0d0,0d0)
        psi2 = (0d0,0d0)
        psi3 = (0d0,0d0)
!        psi4 = (0d0,0d0)
        do j=1,np
        xj = xmin+dx*(j-1)
          psi0 = psi+c(i)*exp(-al/2d0*(x(i)-q(i))**2+
     +             im*p(i)*(x(i)-q(i)))*anrm
          psi1 = psi1+c(i)*exp(-al/2d0*(x(i)+h-q(i))**2+
     +             im*p(i)*(x(i)+h-q(i)))*anrm
          psi2 = psi2+c(i)*exp(-al/2d0*(x(i)+2d0*h-q(i))**2+
     +             im*p(i)*(x(i)+2d0*h-q(i)))*anrm
          psi3 = psi3+c(i)*exp(-al/2d0*(x(i)+3d0*h-q(i))**2+
     +             im*p(i)*(x(i)+3d0*h-q(i)))*anrm

        enddo
        den0 = dlog(abs(psi))
        den1 = dlog(abs(psi1))
        den2 = dlog(abs(psi2))
        den3 = dlog(abs(psi3))
!        den4 = dlog(abs(psi4))
        r0 = (den1-den0)/h
        r1 = (den2-den1)/h
        r2 = (den3-den2)/h
        dr0 = (r1-r0)/h
        dr1 = (r2-r1)/h
        d2r0 = (dr1-dr0)/h
        u(i) = -1d0/4d0/am(i)*(r0**2+dr0)
        du(i) = -1d0/4d0/am(i)*(2d0*r0*dr0+d2r0)
      enddo


      return
      end subroutine
