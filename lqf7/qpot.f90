      subroutine qpot(r0,Nb,ddv,c,am,p,x,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: Nb
      real*8, intent(OUT) :: du(Nb),u(Nb)      
      integer INFO
      real*8 :: w(Nb),x(Nb),c(nb),p(nb),u2(nb),du2(nb),r0(nb,nb)
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
        print *, "QPOT: matrix fails"
        stop
        end if

! the momentum operator r=cf
        do i=1,Ntraj
          rx(i)=cr(1,1)*x(i)+cr(2,1)

! calculate quantum potential
          u(i)=-rx(i)**2/(2d0*am)-cr(1,1)/(2d0*am)
          du(i)=-rx(i)*cr(1,1)/am
        end do

      call rep(nb,w,r0,x,p,u2,du2)

      do i=1,nb
        u(i) = u(i)+u2(i)
        du(i) = du(i)+du2(i)
      enddo

      return
      end subroutine

!------------------------------------------------------
!------extra repulsion added like quantum potential
!----------------------------------------------------
      subroutine rep(nb,w,r0,x,p,u2,du2)
      use cdat
      implicit real*8(a-h,o-z)
      real*8 x(nb),p(nb),dxp(nb),r0(nb,nb),w(nb),u2(nb),du2(nb)
      ww = 0d0
      do i=1,nb
        ww = w(i)+ww
      enddo

      u2 = 0d0

      do i=1,nb
        do j=1,nb
          if(j .ne. i) then
          r = x(i)-x(j)
          dr = r-r0(i,j)
          u2(i) = u2(i)+w(i)/abs(dr)
          du2(i) = du2(i)-dsign(1d0,dr)*w(i)/dr**2
          endif
        enddo 
      enddo

      return
      end subroutine
