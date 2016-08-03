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

      write(*,*) fin
! Open the file
      open(500,file=fin)
     
! Skip the first line 
      read(500,*) 


! Read in the data from the vort files 
      
      do j=1,nz
         do i=1,nx
            read(500,*) xx(i,j),zz(i,j),u(i,j),w(i,j),rhop(i,j),rho(i,j)
         enddo
      enddo   


!      end do 


!*************************************************************
! Calculate integral of the kinetic Engergy in the horizontal 
! for both the u and w velocity fields

! This subroutine is in integral.f 

      call KEh(u,xx,nx,nz,Euu,rho)
      call KEh(w,xx,nx,nz,Eww,rho)
!**************************************************************



!      do i=1,nz
!      write(*,*) Euu(i),Eww(i)
!      enddo 
!      write(*,*) Eww(30)
!      write(*,*) u(30,30), w(30,30), xx(30,30)

!      write(*,*) u(30,30)**2 * (xx(30+1,30)-xx(30-1,30))/2* rho(30,30)

!****************************************************************

! Now that we calculated the integral in the horizontal we will 
! calculate the integral in the vertical, however the points are 
! equally spaced so we will be using a slightly different scheme 


! This scheme is defined in the subroutine integrcol defined in 
! derv.f in the mode1_bbl code 
!***************************************************************

!      do j=1,nz
!         zpts(j)=zz(1,j)
!          write(*,*) zpts(j)
!      enddo 


! we need io build zh because we dont have it any of the files we 
! are including

!      do k=1,nsubd
!      zfrac=1./float(nsubd)
!      zh(k)=zlen*zfrac
!      enddo



C-Set up collocation points      
      call jacobl(nzm,0.,0.,zpts,nzm)

!      do i=0,nzm
!         write(*,*) zpts(i)
!      enddo 

C-Set up weights for use in weak formulation
C-Also set up legendre polynomials
      call quad(nzm,zpts,wg,nzm)



!      do j=1,nzm
!         zpts(j)=zz(1,j)
!          write(*,*) wg(j)
!      enddo
      

      call KEv(Euu,wg,KEu)
      call KEv(Eww,wg,KEw)



!      write(12,*) KEu,KEw


!*************************************************************
! Calculate integral of the Potential Engergy in the horizontal 

! This subroutine is in integral.f 

      call PEh(zz,nx,nz,PEh1,rhop)

!**************************************************************
      call KEv(PEh1,wg,PE)

!      do j=1,nz
!      write(*,*) PEh1(j)
!      enddo

!      call KEv(PEh1,wg,PE)


! Get the time that each file was dumped at 
!      open(69,file='animdump.log')
      read(69,'(37x,F16.16)') time
!      close(69)


      write(12,'(1x,6(F25.16,3x))') time,KEu,Kew,KEu+Kew,PE,
     + KEu+KEw+PE



C-Calculate Integral KEu KEw

!      KEu = 0.0
!      KEw = 0.0

!      do ks=1,nsubd
!        ybar = zh(ks)

!        do kloc = 1,nzloc
!          k = (ks-1)*nzloc + kloc
!          zz = z1(k)

!          a(kloc) = Euu(k)
!          b(kloc) = Eww(k)
c         write(*,*) zz,a(kloc) ! ,fnum,fden
!        enddo

!        outinta = integrcol(a,wg,ybar)
!        outintb = integrcol(b,wg,ybar)
c       write(*,*) nzm,ybar,outint
!        KEu = KEu + outinta
!        KEw = KEw + outintb

!      enddo


!      enddo

      enddo 

! close time file 
      close(69)




      End Program Ecalc




!      Module integration_library 

!      IMPLICIT NONE
      
!      contains
!*********************************************************************
! Integrate the Kinetic Energy in the horizontal diraction

!      subroutine KEh(vel,pos,nx,nz,KEout)
!*********************************************************************

! counters
!      INTEGER :: i,j
! intermediate value
!      real, dimension(nx,nz) :: KEhor
! inputs
!      INTEGER, INTENT(in) :: nx,nz
!      real, dimension(nx,nz), INTENT(IN) :: vel,pos
! output
!      real, dimension(nz), INTENT(out) :: KEout

!*********************************************************************


! Calculate the Integral of velocity velocity squared in the horizontal

!      do j=1,nz
!         do i=1,nx   

! Use forward finite diff for first entry
!            if (i.eq.1) then
!              KEhor(i,j)=1/2*vel(i,j)*vel(i,j)*1/2*(pos(i+1,j)-pos(i,j))

! Use center finite diff for interior values
!            elseif (i.ne.nx) then
!              KEhor(i,j)=1/2*vel(i,j)*vel(i,j)*(pos(i+1,j)-pos(i,j))

! Use backward finite diff for last entry 
!            elseif (i.eq.nx) then
!              KEhor(i,j)=1/2*vel(i,j)*vel(i,j)*1/2*(pos(i+1,j)-pos(i,i))

!            endif
! sum across the horizontal  
!           KEout(j)=KEout(j)+KEhor(i,j)
                    
!         enddo
!      enddo

      
!      end subroutine KEh

!********************************************************************

!      end module integration_library 
      
      
c********************************************************************

      



      
