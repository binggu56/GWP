      program main
      implicit real*8(a-h,o-z)
      real*8, allocatable :: f(:)

      nd = 3
      allocate(f(nd))

      f(1) = 1.0
      f(2) = 2.0
      f(3) = 3.0
      print *,f
      call add(nd,f)
      print *,f
      nd = nd-1
      call add(nd,f)
      print *,f

!      do i=1,nd
!        print *,i
!        if (i==1) i=i+1
!      enddo

!      print *,f(3)
      end program

      subroutine add(nd,f)
      implicit real*8(a-h,o-z)
      real*8 f(nd)
      
      do i=1,nd
        f(i) = f(i)+1d0
      enddo

      return
      end subroutine
