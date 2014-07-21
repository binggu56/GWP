      subroutine qpot(nb,ddv,c,am,p,x,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer,parameter :: ndim = 1
      integer*4,intent(IN) :: nb
      real*8 :: w(nb),x(nb),am,p(nb),ddv(nb),r0(nb),r(nb),dr(nb),ddr(nb)
      complex*16,intent(in) :: c(ntraj)
      real*8, intent(OUT) :: u(ntraj),du(ntraj)
      complex*16 s,z0,z,psi(np),psi0,psi1,psi2,psi3
      real*8 u2(nb),du2(nb),cr(ndim+1)
      common/wave/al,q0,p0,a0
      common/grid/xmin,xmax,dx

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
      h = 2d-4
      do i=1,nb
        psi1 = (0d0,0d0)
        psi0 = (0d0,0d0)
!        psi2 = (0d0,0d0)
!        psi3 = (0d0,0d0)
        do j=1,nb
!          xj = xmin+dx*(j-1)
          psi0 = psi0+c(j)*exp(-al/2d0*(x(i)-x(j))**2+im*p(j)*(x(i)-x(j)))*anrm
          psi1 = psi1+c(j)*exp(-al/2d0*(x(i)+h-x(j))**2+im*p(j)*(x(i)+h-x(j)))*anrm
!          psi2 = psi2+c(j)*exp(-al/2d0*(x(i)+2d0*h-x(j))**2+
!     +             im*p(j)*(x(i)+2d0*h-x(j)))*anrm
!          psi3 = psi3+c(j)*exp(-al/2d0*(x(i)+3d0*h-x(j))**2+
!     +             im*p(j)*(x(i)+3d0*h-x(j)))*anrm
        enddo
!        a0 = dlog(real(conjg(psi0)*psi0))
!        a1 = dlog(real(conjg(psi1)*psi1))
        a0 = dlog(abs(psi0))
        a1 = dlog(abs(psi1))
!        den4 = dlog(abs(psi4))
        r0(i) = (a1-a0)/h
!        r1 = (den2-den1)/h
!        r2 = (den3-den2)/h
!        dr0 = (r1-r0)/h
!        dr1 = (r2-r1)/h
!        d2r0 = (dr1-dr0)/h
!        u(i) = -1d0/4d0/am*(r0**2+dr0)
!        du(i) = -1d0/4d0/am*(2d0*r0*dr0+d2r0)
      enddo

!     fit r with linear basis (1,x)
!      write(*,*) 'x & r', (x(i),i=1,3),(r0(i),i=1,3)
      
      call polyfit(nb,x,r0,ndim,cr)

      do i=1,nb
        r(i) = cr(1)+cr(2)*x(i)
        dr(i) = cr(2)
        ddr(i) = 0d0
        u(i) = -1d0/2d0/am*(r(i)**2+dr(i))
        du(i)= -1d0/2d0/am*(2d0*r(i)*dr(i)+ddr(i))
      enddo
!      do i=1,nb
!        r(i) = cr(1)+cr(2)*x(i)+cr(3)*x(i)**2+cr(4)*x(i)**3
!        dr(i) = cr(2)+2d0*x(i)*cr(3)+3d0*cr(4)*x(i)**2
!        ddr(i) = 2d0*cr(3)+6d0*cr(4)*x(i)
!        u(i) = -1d0/2d0/am*(r(i)**2+dr(i))
!        du(i)= -1d0/2d0/am*(2d0*r(i)*dr(i)+ddr(i))
!      enddo

!      du = 1.1d0*du 

!      call rep(nb,w,r0,x,p,u2,du2)

!     write(*,*) 'du',(du(i),i=1,4)
!      u2 = 0d0
!      du2 = 0d0

!      do i=1,nb
!        u(i) = u(i)+u2(i)
!        du(i) = du(i)+du2(i)
!      enddo

!      write(*,*) (du(i),i=1,4)

      return
      end subroutine



       subroutine polyfit(m,vx,vy, d,c)
!--------------------------------------------------
!      Polynomial fitting of data points 
!      m  : number of points
!      d  : degree of polynomials,up to x**d
!-------------------------------------------------
       implicit real*8(a-h,o-z)
       integer*4, intent(in)       :: m,d
       real*8,dimension(d+1)       :: c 
       integer*4, dimension(d+1)   :: ipiv
       real*8,dimension(d+1)       :: work
       real*8,dimension(m),intent(in) :: vx, vy
    
       real*8 :: X(m,d+1), XT(d+1,m), XTX(d+1,d+1)    
!       integer :: i, j
       integer :: n, lda, lwork
       integer :: info

    
       n = d+1
       lda = n
       lwork = n
    
!       allocate(ipiv(n))
!       allocate(work(lwork))
!       allocate(XT(n, size(vx)))
!       allocate(X(size(vx), n))
!       allocate(XTX(n, n))
    
!     prepare the matrix
       do i = 0, d
          do j = 1, m
             X(j, i+1) = vx(j)**i
          end do
       end do
!      write(*,*) 'polyfit',(vx(i),i=1,3)

    
       XT  = transpose(X)
       XTX = matmul(XT, X)
       c = matmul(XT,vy)

!      print *,x(1,1),x(1,2),x(2,1),x(2,2)
!      print *,'matrix',xtx(1,1),xtx(1,2),xtx(2,1),xtx(2,2)
       ! calls to LAPACK subs DGETRF and DGETRI
       call DGETRF(n, n, XTX, lda, ipiv, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
       call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
       if ( info /= 0 ) then
          print *, 'info =',info," polyfit of r fails"
          return
       end if
       c = matmul( matmul(XTX, XT), vy)

!      call DPOSV('U',n,1,XTX,n,c,n,INFO)    
!       if ( info /= 0 ) then
!          print *,'info =',info, " linear fitting fails."
!          stop
!       end if

      
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
      du2 = 0d0

      do i=1,nb
        do j=1,nb
!          if(j /= i) then
          r = x(j)-x(i)
          if(r .ne. 0d0) then
          r = r/abs(r0(i,j))
!          u2(i) = u2(i)+w(i)/abs(dr)
          du2(i) = du2(i)+dsign(1d0,r)*0.2d0*exp(-4d0*abs(r))/r**2
!          du2(i) = du2(i)-1d0/r**2
          endif
        enddo 
      enddo
!      write(*,*) (du2(i),i=1,6)      
      return
      end subroutine

       
