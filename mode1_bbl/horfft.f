c*********************************************************************
      subroutine horfft (b,is,nk)
c*********************************************************************
      include 'dim.h'
c
      include 'comfft.h'
      dimension b(*),a(nxppy)
c     dimension b(ntot),a(nxppy)
c
      if (is.eq.-1) then
c
c  forward transform, i.e. from physcial to spectral space
c
      do 1000 kz=1,nk
       kk = (kz-1)*nxppy
c
       do i=1,nxppy
        a(i) = b(i+kk)
       enddo
c
         call fft991 (a,work,trigsx,ifaxx,1,nxpp,nx,ny,is)
!-If performing a 3-D simulation will need to do additional
!-complex-to-complex FFT in y-direction
         if (ndims.eq.3) then
           call cfft99 (a,work,trigsy,ifaxy,nxhp,1,ny,nxh,is)
         endif
c
         do 1 i=1,nxppy
            a(i) = scaley*a(i)
    1    continue
c
       do i=1,nxppy
        b(i+kk) = a(i)
       enddo

c      jc = 1
c      ic = (jc-1)*nxpp
c      write(*,*) kz,b(ic+3+kk),b(ic+4+kk)
c
1000   continue
c
      else
c
c  inverse transform, from spectral to physical space
c
      do 2000 kz=1,nk
       kk = (kz-1)*nxppy
c
       do i=1,nxppy
        a(i) = b(i+kk)
       enddo
c
!-If performing a 3-D simulation will need to first do
!-complex-to-complex FFT in y-direction
         if (ndims.eq.3) then
           call cfft99 (a,work,trigsy,ifaxy,nxhp,1,ny,nxh,is)
         endif
         call fft991 (a,work,trigsx,ifaxx,1,nxpp,nx,ny,is)
c
       do i=1,nxppy
        b(i+kk) = a(i)
       enddo
c
2000   continue
c
      endif
c
      return
      end
