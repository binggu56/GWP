       subroutine polyfit(m,vx, vy, d,polyfit)
!--------------------------------------------------
!      Polynomial fitting of data points 
!      m  : number of points
!      d  : degree of polynomials
!-------------------------------------------------
       implicit real*8(a-h,o-z)
       integer, intent(in)         :: d
       dimension(d+1)              :: polyfit,ipiv
       dimension(d+1)              :: work
       dimension(m), intent(in)    :: vx, vy
    
       real*8 :: X(m,d+1), XT(d+1,m), XTX(d+1,d+1)
    
       integer :: i, j
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
       do i = 0, d
          do j = 1, m
             X(j, i+1) = vx(j)**i
          end do
       end do
    
       XT  = transpose(X)
       XTX = matmul(XT, X)
    
       ! calls to LAPACK subs DGETRF and DGETRI
       call DGETRF(n, n, XTX, lda, ipiv, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
       call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
       if ( info /= 0 ) then
          print *, "problem"
          return
       end if
    
       polyfit = matmul( matmul(XTX, XT), vy)
      
       return 
       end subroutine
