      subroutine eigen(N,mat,cnum)
      use cdat
      implicit real*8(a-h,o-z)
      INTEGER,intent(in) :: N
      complex*16, intent(in) :: mat(ntraj,ntraj)
      complex*16 :: A(n,n)
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
      INTEGER          INFO, LWORK
*     .. Local Arrays ..
*     RWORK dimension should be at least MAX(1,3*N-2)
      DOUBLE PRECISION W( N ), RWORK( 3*N-2 )
      COMPLEX*16       WORK( LWMAX )
      LDA = N
*     .. External Subroutines ..
!      EXTERNAL         ZHEEV
!      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
!      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
!      WRITE(*,*)'ZHEEV Example Program Results'
*
*     Query the optimal workspace.
      do j=1,n
        do i=1,n
          A(i,j) = mat(i,j)
        enddo
      enddo
*
      LWORK = -1
      CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF

*     Print eigenvalues.
c      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )
      cnum = w(n)/w(1)

*     Print eigenvectors.
c      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A,
c     $                   LDA )
      
      return
      END subroutine

*     End of ZHEEV Example.
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing a real matrix.
*
      SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,E12.6) )
      RETURN
      END
!--------read matrix from external file--------------
      subroutine read_matrix(nd,A)
      implicit real*8(a-h,o-z)
      integer nd
      complex*16 A(nd,nd)
      open(10,file='matrix')
      do i=1,nd
        read(10,*) (A(i,j),j=1,nd)
      enddo
      close(10)

      return
      end subroutine
