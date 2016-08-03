!*****************************************************
      Program Buildvel
!*****************************************************
! Purpose - build a velocity and density field that 
! has an easily calculated KE and PE 

      use dims
      use input 
!*****************************************************

      IMPLICIT NONE


      INTEGER :: i,j,ii,jj

      real, dimension(:) :: x(nx), z(nz)
      real, dimension(nx,nz) :: xx,zz,u,w,rho 
      
      real :: pi 


      character(21) :: filnam

!*****************************************************
 
      pi = 4 * atan (1.0)


      filnam = 'vortanimout_00001.dat'
 
! read in the gird spacing from the data files

      open(299,file='x_grid.dat')
      open(399,file='z_grid.dat')

      do i = 1,nx
         read(299,*) ii , x(i)
      enddo 

      do j = 1,nz
          read(399,*) jj, z(j)
      enddo 


! build the position fields 

      do j = 1,nz
         do i= 1,nx
            xx(i,j)=x(i)
            zz(i,j)=z(j)
         enddo
      enddo 

! build density field 



      do j = 1,nz
         do i = 1,nx 

            u(i,j) =  sin(x(i)*(2*pi/xlen))
            w(i,j) =  cos(x(i)*(2*pi/xlen))
            rho(i,j) = brunt**2*rho0*z(j)/grav+rho0

         enddo 
      enddo 


      open(200,file=filnam)
      write(200,*) 'ZONE F=POINT, I=',nx,', J=',nz


      do j = 1,nz
         do i = 1,nx
          

         write(200,'(1x,6(F25.16,3x))') xx(i,j),zz(i,j),u(i,j),w(i,j),
     + rho(i,j)

         enddo
      enddo



      end Program Buildvel 
