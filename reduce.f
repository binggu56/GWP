      program main
      implicit real*8(a-h,o-z)
      parameter(nb=5)
      real*8  x(nb),p(nb),c(nb)
 
      x=(/1d0,4d0,5.1d0,4d0,3d0/)
      p = x
      cut = 5d0
      km = nb

10    do i=1,km
        if(x(i) > cut) then
          do j=1,i-1
            x(j) = x(j)
            p(j) = p(j)
          enddo
          km = km-1
          do j=i,km
            x(j) = x(j+1)
            p(j) = p(j+1)
          enddo
          
          goto 10
        endif
      enddo
      
      do i=1,km
        print *,x(i),p(i)
      enddo
      end program

        
