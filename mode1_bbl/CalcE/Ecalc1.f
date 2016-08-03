!- Calculate the phase speed of the wave from the vortanimout_.dat files 

      Program Ecalc
!********************************

        



       use dims
       use integration_library 
          
      IMPLICIT NONE 
     
      INTEGER :: N1,N2
      INTEGER :: t,i,j
!      INTEGER :: nx,ny,nz,nzloc,nsubd,nz2,nzaux,nsdtype,ndims,nxpp
!      INTEGER :: nxh,nxhp,nxf,nyf,nyh,ntot,ntot2,nxppy

      character(12) fspec
      character(9) fmt2
      character(23) fin


!      real :: xlen,zlen
      real, dimension(nx,nz) :: xx,zz,u,w,rhop,rho
!      real, allocatable, dimension(:,:) :: Eu,Ew
      real, dimension(nz) :: Euu,Eww
!      real, dimension(nsubd) :: zh
!      real :: ybar

!      real, dimension(nsubd*nzloc) :: a,b


      write(*,*) "Enter the first and last vort files to be used"
      read(*,*) N1,N2
      

      fmt2 = '(i5.5,a4)'
      fspec = 'vortanimout_'



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
      
      do j=1,nx
         do i=1,nz
            read(500,*) xx(i,j),zz(i,j),u(i,j),w(i,j),rhop(i,j),rho(i,j)
         enddo
      enddo   


      end do 

      call KEh(u,xx,nx,nz,Euu,rho)
      call KEh(w,xx,nx,nz,Eww,rho)


      do i=1,nz
      write(*,*) Euu(i),Eww(i)
      enddo 
!      write(*,*) Eww(30)
!      write(*,*) u(30,30), w(30,30), xx(30,30)

!      write(*,*) u(30,30)**2 * (xx(30+1,30)-xx(30-1,30))/2* rho(30,30)

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

      



      
