      subroutine reduce(nb,x,p,c)
      implicit real*8(a-h,o-z)
      real*8  x(nb),p(nb),c(nb)
      real*8 ax(nb),ap(nb),ac(nb)
      cut = 5d0
      ax = 0d0
      ap = 0d0
      ac = 0d0

10    do i=1,nb
        if(x(i) > cut) then
          do j=1,i-1
            x(j) = x(j)
            p(j) = p(j)
          enddo
          nb = nb-1
          do j=i,nb
            x(j) = x(j+1)
            p(j) = p(j+1)
          enddo
          
          goto 10
        endif
      enddo



        
