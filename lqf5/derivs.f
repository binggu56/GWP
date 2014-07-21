!----------potential---------------------
      subroutine derivs(x,Nb,v,dv,ddv)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Nb
      real*8,intent(IN) :: x(Ntraj)
      real*8,intent(OUT) :: v(Ntraj),dv(Ntraj),ddv(Ntraj)

      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx

!------harmonic----------------------
!      ak=1.d0
!      do i=1,Ntraj
!        dv(i) = ak*x(i)
!        v(i)  = ak*x(i)**2/2d0
!        ddv(i) = ak
!      enddo
!------morse potential---------------
      de = 0.176d0
      x0 = 1.4d0
      a = 1.02d0
!      do i=1,Nb
!        d = (1d0-exp(-a*(x(i)-x0)))
!        v(i) = de*d**2
!        dv(i) = 2*de*d*a*exp(-a*(xj-x0))
!        ddv(i)=2*de*(-d*exp(-a*((x(i)-x0)))*a**2+
!     +         (exp(-a*(x(i)-x0)))**2*a**2)
!      enddo

      av = 8d0
      ax = 2d0
      do i=1,nb
        v(i) = av*exp(-ax*x(i)**2)
        dv(i) = -2d0*ax*av*x(i)*exp(-ax*x(i)**2)
        ddv(i) = av*(-2d0*ax*exp(-ax*x(i)**2)+4d0*x(i)**2*
     +           ax**2*exp(-ax*x(i)**2))
      enddo

!---------doulbe well---------------------------
!      eta = 1.3544d0
!      do i=1,Ntraj
!        v(i) = 1d0/16d0/eta*x(i)**4-0.5d0*x(i)**2
!        dv(i) = 1d0/4d0/eta*x(i)**3-x(i)
!        ddv(i) = 3d0/4d0/eta*x(i)**2-1d0
!      enddo

      return
      end subroutine

