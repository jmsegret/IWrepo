
      module integration_library

       use input
!      use comleng


      IMPLICIT NONE

      contains
!*********************************************************************
! Integrate the Kinetic Energy in the horizontal diraction

      subroutine KEh(vel,pos,nx1,nz1,KEout,rho1)
!*********************************************************************


! inputs
      INTEGER, INTENT(in) :: nx1,nz1
      real, dimension(nx1,nz1), INTENT(IN) :: vel,pos,rho1
! output
      real, dimension(nz1), INTENT(out) :: KEout

! counters
      INTEGER :: i,j
! intermediate value
      real, dimension(nx1,nz1) :: KEhor
      real :: dp , rv , KEe
!*********************************************************************


! Calculate the Integral of velocity velocity squared in the horizontal
!          j=1
       do j=1,nz1
         do i=1,nx1

! Use forward finite diff for first entry
!            if (i.eq.1) then
!            dp = (pos(i+1,j)-pos(i,j))/2
!            rv =rho1(i,j)*vel(i,j)**2
!            KEhor(i,j) =  rv!*dp
!             KEe=rho1(i,j)*vel(i,j)**2


!              KEhor(i,j)=1/2*vel(i,j)*vel(i,j)*1/2*(pos(i+1,j)-pos(i,j))

! Use center finite diff for interior values
!            else!if (i.lt.nx1) then
!            dp = (pos(i+1,j)-pos(i-1,j))/2
!            rv =(rho1(i,j)+rho0-1.)/2.*vel(i,j)**2
             rv=rho0/2.*vel(i,j)**2

            KEhor(i,j) =  rv!*dp   

    


!        KEhor(i,j)=1/2*(vel(i,j)**2)*(pos(i+1,j)-pos(i-1,j))/2*rho1(i,j)
!              KEhor(i,j)=vel(i,j)
!               write(*,*) KEhor(i,j),dp,rv 
!              write(*,*) KEhor(i,j),vel(i,j),(pos(i+1,j)-pos(i-1,j))/2
! Use backward finite diff for last entry 
!            elseif (i.eq.nx1) then
!            dp = (pos(i,j)-pos(i-1,j))/2
!            rv =rho1(i,j)*vel(i,j)**2
!            KEhor(i,j) =  rv*dp


!             KEe=(KEe+rv)


!              KEhor(i,j)=1/2*vel(i,j)*vel(i,j)*1/2*(pos(i+1,j)-pos(i,i))

!            endif
! sum across the horizontal  
!           KEout(j)=KEout(j)+KEhor(i,j)

         enddo
!         write(*,*) KEhor(i,j)         
         KEout(j) =(xlen)*sum (KEhor(:,j))/nx1
      enddo
       
!       write(*,*) KEout
      end subroutine KEh

!*******************************************************************

      Subroutine KEv(KEh,wg1,KEout)

!*******************************************************************

! Purpose - calculate the integral of the kinetic energy in the vertical 

!************************************************************************
       use dims
       use input, only: zlen

! inputs
!      INTEGER, INTENT(in) :: nz1,nzm,nsubd
      real, dimension(1:nz), INTENT(IN) :: KEh
      real, dimension(nzm), INTENT(IN) :: wg1
      real, dimension(nsubd) :: zh
! output
      real, INTENT(out) :: KEout

! counters
      INTEGER :: k,ks,kloc
! intermediate value
      real :: KE,ybar,zfrac,zfrac1,zn 
      real, dimension(nzloc) :: a
!*********************************************************************

!      real, dimension(*) :: amat1(1:nzloc),wg1(0:nzm)
!      real :: sum,ybar1
!      Integer :: k
      open(19,file='subd.dat')
      read(19,*)

      do k=1,nsubd

      read(19,*) zn, zfrac1

      zfrac=zfrac1/zlen

!      zfrac=1./float(nsubd)
      zh(k)=zlen*zfrac

!      write(*,*) 'zh =', zh(k)
      enddo
      
      close(19)      

      KEout = 0.0 


      do ks=1,nsubd
        ybar = zh(ks)

        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          a(kloc) = KEh(k)
!          write(*,*) 'KEh = ' , KEh(k)
        enddo



      KE = 0.0

      do k=1,nzloc 

!          write(*,*) wg1(k)
!          KE = KE + wg1(k)*a(k)
!this is what was here i am trying to find a bug 
          KE = KE + wg1(k-1)*a(k)
      enddo
 

c     write(*,*) sum,ybar
      KE = KE*ybar/2.
       KEout = KEout + KE

       enddo 

!       KEout = KEout + KE

      end subroutine KEv







!*********************************************************************
! Integrate the Availabe Potential  Energy in the horizontal diraction

      subroutine PEh(pos,nx1,nz1,KEout,rho1)
!*********************************************************************


! inputs
      INTEGER, INTENT(in) :: nx1,nz1
      real, dimension(nx1,nz1), INTENT(IN) :: pos,rho1
! output
      real, dimension(nz1), INTENT(out) :: KEout

! counters
      INTEGER :: i,j
! intermediate value
      real, dimension(nx1,nz1) :: PEhor
      real, dimension(nz1) :: rhobar
      real :: dp , pe , KEe
! parameters
!      real, parameter :: g=9.81 
!*********************************************************************



! Calculate the Potential Enegergy of the field. 

       do j=1,nz1
         do i=1,nx1

! Calculate the rhobar(z) from the brunt
! brunt and rho0 are both from input.f 

       rhobar(j)= -(brunt**2)*rho0*pos(i,j)/grav+1.



      
! Subtract out the background density 
! So we can calculate the availabe PE 

!            pe =grav*(rho1(i,j)-rhobar(j))*pos(i,j)   
           pe=1./2./rho0*grav**2/brunt**2*rho1(i,j)**2!*pos(i,j) 
           PEhor(i,j) =  pe

         enddo
!         write(*,*) 'z =', pos(i,j), 'rhop = ',rho1(i,j)
         KEout(j) =(xlen)*(sum(PEhor(:,j))+PEhor(1,j))
     +             /((nx1))
      enddo

!       write(*,*) KEout
      end subroutine PEh

!***********************************************************************


!*******************************************************************
























!********************************************************************

      end module integration_library


c********************************************************************

