       subroutine testderv
C-Tests derivative evaluation scheme
C-2/14/03: In this case the conventional
C-sequential differentiation concept is tested.
       logical cbch

       include 'parambctype.h'
       include 'dim.h'

       include 'comcoeff.h'
       include 'comleng.h'
       include 'comwave.h'
       include 'comsubd.h'
       include 'comgridptxyz.h'

       parameter (tiny = 1.e-30)   

       dimension f(nz),df(nz),d2f(nz),d3f(nz)
       dimension tmp1(nzloc),tmp2(nzloc),tmp3(nzloc),
     >           tmp4(nzloc)

       write(*,*) 'Testing derivative calculation'

       open(700,file='testderv.dat',status='unknown')
       write(700,*) 'VARIABLES = "z","f","Df Num","D^2f Num","D^3f Num",
     >                           "Df Anal","D^2f Anal","D^3f Anal",
     >                           "Error in D","Error in D^2",
     >                           "Error in D^3"'

       open(710,file='testin.dat',status='old')
       rewind(710)

       do 10 ks=1,nsubd
c      do 10 ks=nsubd,1,-1
C-Set Up function
         ybar = zh(ks)
c        write(*,*) 'YBAR',zh(ks)
         do kloc = 1,nzloc
c        do kloc = nzloc,1,-1
           k = (ks-1)*nzloc + kloc
c          f(k) = fin(z(k))
           read(710,*) dummy,f(k),dummy
           write(*,*) f(k)
           tmp1(kloc) = f(k) ! fin(z(k))
         enddo

C-Calculate 1st derivative
         call dpdzcol(tmp1,tmp2,d,zpts,nzm,ybar,1,nzloc)     
C-Calculate 2nd derivative
c        call dpdzcol(tmp2,tmp3,d,zpts,nzm,ybar,1,nzloc)
         call d2pdz2col(tmp1,tmp3,d2,zpts,nzm,ybar,1,nzloc)
C-Calculate 3rd derivative
c        call dpdzcol(tmp3,tmp4,d,zpts,nzm,ybar,1,nzloc)
         call d3pdz3col(tmp1,tmp4,d3,zpts,nzm,ybar,1,nzloc)

C-Update functions
         do kloc = 1,nzloc
           k = (ks-1)*nzloc + kloc
c          write(*,*) z(k),tmp2(kloc) ! ,tmp2(kloc),tmp3(kloc),tmp4(kloc)
           df(k) = tmp2(kloc)
           d2f(k) = tmp3(kloc)
           d3f(k) = tmp4(kloc)
         enddo
 10    continue

C-Output value and calculate relative errors
       do k=1,nz
         errord1 = abs( (df(k) - danal(z(k))) /
     >                  (danal(z(k))+tiny))
         errord2 = abs( (d2f(k) - d2anal(z(k)))/
     >                  (d2anal(z(k))+tiny))
         errord3 = abs( (d3f(k) - d3anal(z(k)))/
     >                  (d3anal(z(k))+tiny))
     
c        if (danal(z(k)).eq.0.) errord1 = 0.
c        if (d2anal(z(k)).eq.0.) errord2 = 0.
c        if (d3anal(z(k)).eq.0.) errord3 = 0. 

         write(700,'(1x,F10.5,10(2x,E12.6)))')
     >                z(k),f(k),df(k),d2f(k),d3f(k),
     >                danal(z(k)),d2anal(z(k)),d3anal(z(k)),
     >                errord1,errord2,errord3
       enddo 
      
       close(700)

       return
       end


C************************************************************

       function fin(z)
       z00 = 0.5
c      fin = a*exp(-(z-0.5)**2./b)
       rr2=2.772588722*((z-z00)**2)
       r3=0.2
       umax = 3.2
       fin = umax*0.375*exp(-rr2/r3**2.)
c      fin = (z-0.5)**2.
       write(*,*) fin
       return
       end

       function danal(z)
       a=1e10 
       b=0.005
       danal = -2.*(z-0.5)*a*exp(-(z-0.5)**2./b)/b
       return
       end

       function d2anal(z)
       pi = 4.0*atan(1.)
       a = 1. 
       om = 2.*pi*a
       d2anal = -(om**2.)*cos(om*z)
       return
       end

       function d3anal(z)
       pi = 4.0*atan(1.)
       a = 1. 
       om = 2.*pi*a
       d3anal = (om**3.)*sin(om*z)
       return
       end 
   
 

