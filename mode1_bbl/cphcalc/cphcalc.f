!- Calculate the phase speed of the wave from the vortanimout_.dat files 

      Program cphcalc
!********************************


      IMPLICIT NONE

      
!      INTEGER :: nsubd      !number of subdomains
!      REAL :: Htotal     !height of domain
!      REAL :: xx,zz,u,w,rhop,rho
      INTEGER :: I , J , N1,N2 , skip , ianimc , x,t, z1,time
      character*12 fspec
      character*9 fmt2
      character*15 fout
      character*23 fin
      character*80 header

      real :: zz,w,rhop,rho,xlen
      real, allocatable,  dimension(:) :: xx,u,derv,time1
      real, allocatable, dimension(:) :: zero,x1,x2,xx1,xx2

      write(*,*) "Input number of vertical and horizontal points as x,y"
      read(*,*) I,J
      write(*,*) "Enter the first and last vort files to be used"
      read(*,*) N1,N2
      write(*,*) "Enter the row of data we want and, horizontal length"
      read(*,*) skip,xlen

      allocate(xx(I),u(I),derv(I),time1(I))
      allocate(zero(2))
      allocate(x1(N2-N1),x2(N2-N1),xx1(N2-N1),xx2(N2-N1))

      open(399,file='animdump.log')

      if (N1.ne.1) then
         do t=1,N1-1
            read(399,*)
         enddo 
      endif 

      do t=N1,N2
         read(399,'(37X,1F20.15)') time1(t-N1+1)
         
      enddo 


      open(99,file='zeros.dat')
! We are going from the first file specified 
! to the last file 


      do t=N1,N2

      ianimc = t
      fmt2 = '(i5.5,a4)'

      fspec = 'vortanimout_'
      fin = fspec
      fout = 'vec_'

! Write the name of the file to be opened and the name of the 
! file to be written to 

      write(unit=fin(13:23),fmt=fmt2) ianimc,'.dat'
      write(unit=fout(5:15),fmt=fmt2) ianimc,'.dat'

      open(500,file=fin)
      open(200,file=fout)
     

      read(500,*) 


! This is the header needed for tecplot 

      write(200,*) 'ZONE F=POINT, I=',I,', J=',1


! Here we are skipping all of the rows in the vort files
! until we reach the first row that corresponds to the row of 
! u velocity data that we are interested in
  
      if (skip.ne.1) then 
      do x=1,I*(skip-1)
         read(500,*)
      enddo
      
 
      endif 
    
! Read in the desired information from the vort files
! store the x location data and u velocity data into arrays      
 
      do x=1,I
         read(500,*) xx(x),zz,u(x),w,rhop,rho
         
      enddo  
! Calculate the derivative of velocity with respect to time
! So that we can use the zero points to find the x location of 
! the maximum and minimum values of the u velocity field 
      do x=1,I   
         if (x.eq.1) then
            derv(x)=(u(x+1)-u(x))/(xx(x+1)-xx(x))
         elseif (x.ne.I) then
            derv(x)=(u(x+1)-u(x-1))/(xx(x+1)-xx(x-1))
         elseif (x.eq.I) then
            derv(x)=(u(x)-u(x-1))/(xx(x)-xx(x-1))
         endif 
            
! write the data to the vec_ file in a format usable by tecplot 
         write(200,'(1x,4(e14.7,3x))') xx(x),zz,u(x),derv(x)
         
      enddo
      
      
! open a new file and call it zeros.dat, this is where we will calculate 
! the location of the zero points of our derivative data 

! First we fount the location of the zero crossing by checking each
! set of neighboring locations to see if the zeros crossing occured there

! then the linear interpolation equation between two points was used
! to find a more accurate location fo the zero crossing 



! The first zero crossing found is put in the first collumn 
! the second one found is put in the second collumn. 

! This is a problem because we are then tracking both zero 
! crossings in both collumns, (Not keeping the location of the
! positive and negative lobes seperate) 
      z1=1

      do x=1,I
         if (derv(x)*derv(x+1).lt.0) then

             
            zero(z1)=derv(x)*(xx(x+1)-xx(x))/(derv(x+1)-derv(x))+xx(x)
            z1=z1+1
            
         endif
      enddo 
       
! A time is needed to find the phase speed.
! Here a 5 second offset between each vec file is used
! this is because that is the time for dumping our a vort file
! Although I am not sure how accurate the dumping is ( if the 
! code actually dumps out exactly every 5 seconds) 

      write(99,*) time1(N1-t-1), zero
      
!     write(99,file='zero.dat') zero

      enddo

      close(99)
! I close and then open the zeros.dat file so that I 
! can start reading the file from the top 
      open(99,file='zeros.dat')
      open(299,file='cph.dat')
      



 

! Calculate phase speed and write it to the cph.dat file 
! for the combined lobe phase speed 

      do x=1,N2-N1
         read(99,*)  time, x1(x) , x2(x)
      enddo 
 


      

      do x=1,N2-N1+1
         if(x.ne.1) then
            if (x1(x).gt.x1(x-1)) then
!        write(299,*)time1(x),(x1(x+1)-x1(x-1))/(time1(x+1)-time1(x-1)),
!     + (x2(x+1)-x2(x-1))/(time1(x+1)-time1(x-1)),x1(x),x2(x)
            elseif (x1(x).lt.x1(x-1)) then
! flip the two columns from the location where the lobe on the right 
! moves into the field of view on the left side

! Do this in order to keep track of the two lobes seperately
! Also add on the size of the domain so it as if we are tracking the lobe
! constantly moving to the right   
               xx1=x1
               xx2=x2
               x1(x:N2-N1)=xx2(x:N2-N1)
               x2(x:N2-N1)=xx1(x:N2-N1)+xlen
!      write(299,*)time1(x),(x1(x+1)-xx1(x-1))/(time1(x+1)-time1(x-1)),
!     + (x2(x+1)-xx2(x-1))/(time1(x+1)-time1(x-1)),x1(x),x2(x)
             endif 
         endif
      enddo 

      do x=2,N2-N1-1
         write(299,'(1x,5(F25.16,3x))')time1(x),
     + (x1(x+1)-x1(x-1))/(time1(x+1)-time1(x-1)),
     + (x2(x+1)-x2(x-1))/(time1(x+1)-time1(x-1)),x1(x),x2(x)
      enddo

          

! Deallocate the variables that were allocated
! so that we free the memory for later use

      deallocate(xx,u,derv,zero,x1,x2,time1)


      End Program cphcalc
 
      
