      program plynomalFitting
      implicit real*8(a-h, o-z)
     
      ! let us test it
      integer, parameter      :: degree = 2
      integer                 :: i
      real*8, dimension(11) :: x = (/ (i,i=0,10) /)
      real*8, dimension(11) :: y = (/ 1,   6,  17,  34, 
     +                                  57,  86, 121, 162, 
     +                                  209, 262, 321 /)
      real*8, dimension(degree+1) :: a
     
      call polyfit(11, x, y, degree,a)
     
      write (*, '(F9.4)') a
     
      end program

  
