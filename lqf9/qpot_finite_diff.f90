      subroutine qpot(nb,ddv,c,am,p,x,w,u,du)
      use cdat
      implicit real*8(a-h,o-z)
      integer*4,intent(IN) :: nb
      real*8 :: w(nb),x(nb),am,p(nb),ddv(nb),r0(nb),r(nb),dr(nb),ddr(nb)
      complex*16,intent(in) :: c(ntraj)
      real*8, intent(OUT) :: u(ntraj),du(ntraj)
      complex*16 s,z0,z,psi(np),psi0,psi1,psi2,psi3
      real*8 cr(2)
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
      h = 1d-3
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
        den0 = dlog(real(conjg(psi0)*psi0))
        den1 = dlog(real(conjg(psi1)*psi1))
!        den2 = dlog(abs(psi2))
!        den3 = dlog(abs(psi3))
!        den4 = dlog(abs(psi4))
        r0(i) = (den1-den0)/2d0/h
!        r1 = (den2-den1)/h
!        r2 = (den3-den2)/h
!        dr0 = (r1-r0)/h
!        dr1 = (r2-r1)/h
!        d2r0 = (dr1-dr0)/h
!        u(i) = -1d0/4d0/am*(r0**2+dr0)
!        du(i) = -1d0/4d0/am*(2d0*r0*dr0+d2r0)
      enddo

!     fit r with linear basis (1,x)
      write(*,*) (x(i),i=1,3),(r0(i),i=1,3)
      write(*,*) 'nb =',nb
      
      call polyfit(nb,x,r0,1,cr)

      do i=1,nb
        r(i) = cr(1)+cr(2)*x(i)
        dr(i) = cr(2)
        ddr(i) = 0d0
        u(i) = -1d0/2d0/am*(r(i)**2+dr(i))
        du(i)= -1d0/2d0/am*(2d0*r(i)*dr(i)+ddr(i))
      enddo
        
      return
      end subroutine



       subroutine polyfit(m,vx,vy, d,c)
!--------------------------------------------------
!      Polynomial fitting of data points 
!      m  : number of points
!      d  : degree of polynomials
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
    
       ! prepare the matrix
!       do i = 0, d
!          do j = 1, m
!             X(j, i+1) = vx(j)**i
!          end do
!       end do
      write(*,*) (vx(i),i=1,3)
      print *,'m=',m

       do i=1,m
         x(i,1) = 1d0
         x(i,2) = vx(i)
       enddo
         
    
!       XT  = transpose(X)
!       XTX = matmul(XT, X)
!       c = matmul(XT,vy)

      print *,x(1,1),x(1,2),x(2,1),x(2,2)
      print *,xtx(1,1),xtx(1,2),xtx(2,1),xtx(2,2)
       ! calls to LAPACK subs DGETRF and DGETRI
       call DGETRF(n, n, XTX, lda, ipiv, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
       call DGETRI(n, XTX, lda, ipiv, work, lwork, info)

!      call DPOSV('U',2,1,XTX,2,c,2,INFO)    
       if ( info /= 0 ) then
          print *,'info =',info, " linear fitting fails."
          stop
       end if
       c = matmul( matmul(XTX, XT), vy)
      
       return 
       end subroutine



       
