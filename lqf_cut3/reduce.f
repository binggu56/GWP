      subroutine reduce(nb,w,x,p,c)
      use cdat
      implicit real*8(a-h,o-z)
      real*8  x(ntraj),p(ntraj),w(ntraj)
      complex*16 c(ntraj)

      cut = 2.5d0

10    do i=1,nb
        if(x(i) > cut) then
          do j=1,i-1
            x(j) = x(j)
            p(j) = p(j)
            c(j) = c(j)
            w(j) = w(j)
          enddo
          nb = nb-1
          do j=i,nb
            x(j) = x(j+1)
            p(j) = p(j+1)
            c(j) = c(j+1)
            w(j) = w(j+1)
          enddo
          
          goto 10
        endif
      enddo

      return 
      end subroutine

        
