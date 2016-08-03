      module integral 

      use dims, only: nzloc,nzm
      use comleng, only: wg

      IMPLICIT NONE

      contains
c-k: Local GLL grid-point
      real function integrcol(amat1,wg1,ybar1) 

      real, dimension(*) :: amat1(1:nzloc),wg1(0:nzm)
      real :: sum,ybar1
      Integer :: k

      

      sum = 0.0

      do k=1,nzloc
          sum = sum + wg1(k-1)*amat1(k)
      enddo

c     write(*,*) sum,ybar
      integrcol = sum*ybar1/2.
      return
      

      end function
      end module integral

