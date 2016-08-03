!- Calculate the phase speed of the wave from the vortanimout_.dat files 

      Program Ecalc
!********************************

        



       include 'dim.h'
       include 'comleng.h'
       include 'comgridptxyz.h'



      

      
!      INTEGER :: nsubd      !number of subdomains
!      REAL :: Htotal     !height of domain
!      REAL :: xx,zz,u,w,rhop,rho
      INTEGER :: II,JJ,N1,N2,skip,ianimc,i,t,j
      INTEGER :: ks,kloc,k 
      INTEGER :: nx,ny,nz,nzloc,nsubd,nz2,nzaux,nsdtype,ndims,nxpp
      INTEGER :: nxh,nxhp,nxf,nyf,nyh,ntot,ntot2,nxppy

      character*12 fspec
      character*9 fmt2
      character*15 fout
      character*23 fin
      character*80 header

      real :: xlen,zlen
      real, allocatable,  dimension(:,:) :: xx,zz,u,w,rhop,rho
      real, allocatable, dimension(:,:) :: Eu,Ew
      real, allocatable, dimension(:) ::z1,Euu,Eww
      real, dimension(nsubd) :: zh,ybar

      write(*,*) "Input number of vertical and horizontal points as x,z"
      read(*,*) II,JJ
      write(*,*) "Enter the first and last vort files to be used"
      read(*,*) N1,N2
      write(*,*) "Enter the row of data we want and, xlen,zlen"
      read(*,*) skip,xlen,zlen

      allocate(xx(II,JJ),zz(II,JJ),u(II,JJ))
      allocate(w(II,JJ),rhop(II,JJ),rho(II,JJ))
      allocate(Eu(II,JJ),Ew(II,JJ))
      
      allocate(z1(nsubd*nzloc),Euu(JJ),Eww(JJ))

!      open(399,file='animdump.log')

!      if (N1.ne.1) then
!         do t=1,N1-1
!            read(399,*)
!         enddo 
!      endif 

!      do t=N1,N2
!         read(399,'(37X,1F20.15)') time1(t-N1+1)
         
!      enddo 


! We are going from the first file specified 
! to the last file 


      do t=N1,N2

      ianimc = t
      fmt2 = '(i5.5,a4)'

      fspec = 'vortanimout_'
      fin = fspec
!      fout = 'vec_'

! Write the name of the file to be opened and the name of the 
! file to be written to 

      write(unit=fin(13:23),fmt=fmt2) ianimc,'.dat'
!      write(unit=fout(5:15),fmt=fmt2) ianimc,'.dat'

      open(500,file=fin)
!      open(200,file=fout)
     

      read(500,*) 


! This is the header needed for tecplot 

!      write(200,*) 'ZONE F=POINT, I=',I,', J=',J
    
! Read in the desired information from the vort files
! store the data into arrays
      
      do j=1,JJ
         do i=1,II
            read(500,*) xx(i,j),zz(i,j),u(i,j),w(i,j),rhop(i,j),rho(i,j)
         enddo
      enddo   

! Calculate the Integral of velocity velocity squared in the horizontal

      do j=1,JJ
         do i=1,II   
            if (i.eq.1) then
               Eu(i,j)=1/2*u(i,j)*u(i,j)*1/2*(xx(i+1,j)-xx(i,j))
               Ew(i,j)=1/2*w(i,j)*w(i,j)*1/2*(xx(i+1,j)-xx(i,j))
            elseif (i.ne.I) then
               Eu(i,j)=1/2*u(i,j)*u(i,j)*(xx(i+1,j)-xx(i,j))
               Ew(i,j)=1/2*w(i,j)*w(i,j)*(xx(i+1,j)-xx(u,j))
            elseif (i.eq.I) then
               Eu(i,j)=1/2*u(i,j)*u(i,j)*1/2*(xx(i+1,j)-xx(i,i))
               Ew(i,j)=1/2*w(i,j)*w(i,j)*1/2*(xx(i+1,j)-xx(i,j))
            endif 
           Euu(j)=Euu(j)+Eu(i,j)
           Eww(j)=Eww(j)+Ew(i,j)           
         enddo
      enddo

      
! we need to build zh because we dont have it any of the files we 
! are including


      do k=1,nsubd
      zfrac=1./float(nsubd)
      zh(k)=zlen*zfrac
      enddo



! Now we need to integrate in the vertical direction

      


C-Calculate Integral KEu KEw

      KEu = 0.0
      KEw = 0.0

      do ks=1,nsubd
        ybar = zh(ks)

        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          zz = z1(k)

          a(kloc) = Euu(k)
          b(kloc) = Eww(k)
c         write(*,*) zz,a(kloc) ! ,fnum,fden
        enddo

        call integrcol(a,outinta,nzm,ybar)
        call integrcol(b,outintb,nzm,ybar)
c       write(*,*) nzm,ybar,outint
        KEu = KEu + outinta
        KEw = KEw + outintb

      enddo


      enddo
! Deallocate the variables that were allocated
! so that we free the memory for later use

      deallocate(xx,zz,u,w,rhop,rho,z1)



      End Program Ecalc


      
      
c********************************************************************

c-Using properties of Legendre interpolants calculates 
c-integral of a given function in interval of length ybar.
C-ONLY WITHIN A GIVE VERTICAL SUBDOMAIN
C-(Created on 4/2/04 by PD)
c
c
c-a: The input column of data
c-b: The output value for integral
c-nzbot: Bottom point of point range
c-nztop: Top point of point range
      include 'dim.h'
      include 'comleng.h'
      include 'comgridptxyz.h'

      real, dimension :: a(*)
c
C-Set up values of collocation points

c
c-k: Local GLL grid-point
      sum = 0.0

      do k=1,nzloc
          sum = sum + wg(k-1)*a(k)
      enddo

c     write(*,*) sum,ybar
      outint = sum*ybar/2.
c
      return
      end
 


      



      
