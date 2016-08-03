!- Calculate the phase speed of the wave from the vortanimout_.dat files 

      Program Ecalc
!********************************

        



       use dims
       use integration_library 
       use comleng 
  
      IMPLICIT NONE 
     
      INTEGER :: N1,N2
      INTEGER :: t,i,j
!      INTEGER :: nx,ny,nz,nzloc,nsubd,nz2,nzaux,nsdtype,ndims,nxpp
!      INTEGER :: nxh,nxhp,nxf,nyf,nyh,ntot,ntot2,nxppy

      character(12) fspec
      character(9) fmt2
      character(23) fin
      character(21) outfil


!      real :: xlen,zlen
      real, dimension(nx,nz) :: xx,zz,u,w,rhop,rho
      real :: KEu,KEw,PE,time
      real, dimension(nz) :: Euu,Eww,PEh1

      real, dimension(nzm) :: wg
      real, dimension(nz) :: zpts
!      real, dimension(nsubd) :: zh
!      real :: ybar

!      real, dimension(nsubd*nzloc) :: a,b

       pi2 = 2

       xomega = brunt*xkcrit/sqrt(xkcrit**2+xmcrit**2)

       twave = pi2/xomega



! Open file that contains dump times 
      open(69,file='animdump.log')


      write(*,*) "Enter the first and last vort files to be used"
      read(*,*) N1,N2
      
      write(*,*) nx,nz      


      fmt2 = '(i5.5,a4)'
      fspec = 'vortanimout_'
      outfil = 'outfile_'
      write(unit=outfil(8:21),fmt='(i4.4,a1,i4.4,a4)') N1,'_',N2,'.dat'
      open(12,file=outfil)

      do t=N1,N2

      fin = fspec

! Write the name of the file to be opened and the name of the 
! file to be written to 

      write(unit=fin(13:23),fmt=fmt2) t,'.dat'


      End Program Ecalc

      



      
