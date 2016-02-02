       subroutine testhelm
C-Tests generic Helmholtz solver
       logical cbch

       parameter (tiny=1.e-30)
       include 'parambctype.h'
       include 'dim.h'

       include 'comcoeff.h'
       include 'comleng.h'
       include 'comwave.h'
       include 'comsubd.h'
       include 'comgridptxyz.h'
c      include 'comparampen.h'
 
       dimension rhsinr(nz),rhsini(nz)
       dimension ar(nz),error(nz) 
       dimension fr(nz),fi(nz)
       dimension codesol(nz)

c      complex ucol(nz),vcol(nz),wcol(nz),divvel(nz)
c      dimension tr(nz),ti(nz),tlr(nz),tli(nz) 

       open (705,file = 'testlaplace.dat')
       open (700,file='testout.dat',status='unknown')
       write(700,*) 'VARIABLES = "Z","FR","ANALYTIC","ERROR",
     > "|ERROR|"'
!      open (710,file='testin.dat',status='old')
!      rewind(710)
c      read (710,*) 
       pi = 4.0*atan(1.)

         
       write(*,*) 'CALLING TESTHELM'

c      do k=1,nz
c               divvel(k) = cmplx(0.,0.)
c               ucol(k) = cmplx(0.,0.)
c               vcol(k) = cmplx(0.,0.)
c               wcol(k) = cmplx(z(k)*z(k)*0.5,0.)
c      enddo

c      xwin = 0.    
c      ywin = 0.   
c      call divcolc(ucol,vcol,wcol,xwin,ywin,divvel)

c      do k=1,nz
c      write(*,*) z(k),real(wcol(k)),real(divvel(k))
c      enddo

c      pause

C-Set up values of collocation points
C      do k=1,nz
C        z(k) = 0.5*ybar*(zpts(k-1)+1.0)
C      enddo

C-Set up beta coefficient
       coeffbeta =  1.e2
       epsilon = 1./coeffbeta
       write(*,*) 'epsilon',epsilon
       call penalty_setup(epsilon)
 
C-Set up RHS (imaginary and real)
C-for each point in element. This is done using
C-GLOBAL indexing system
       do k=1,nz
         zz = z(k)
         rhsinr(k) = -1.
         rhsini(k) = -1.
       enddo

C-Average at interfaces
c      do ks=2,nsubd
c        k=(ks-1)*nzloc+1
c        rhsavg = 0.5*(rhsinr(k)+rhsinr(k-1))
c        rhsinr(k) = rhsavg
c        rhsinr(k-1) = rhsavg
c      enddo
         
C-Set Up GLOBAL (!!) Boundary Conditions and types
       ifbcbot = idirichlet
!      ifbcbot = ineumann
       ifbctop = idirichlet
!      ifbctop = ineumann
 
       bcrbot = 0.
       bcibot = 0.
       bcrtop = 0.
       bcitop = 0.


       ifuncreal=0
C-Call penalty method for calculation of singularly perturbed
C-boundary value problem (Method borrowed from Hesthaven, 97)
       call penalty_main(epsilon,ifuncreal,
     >                   bcrbot,bcibot,bcrtop,bcitop,
     >                   ifbcbot,ifbctop,rhsinr,rhsini,
     >                   fr,fi)


C-Check what penalty method yields
c      write(*,*) 'Checking Penalty Method Results'
c      call penalty_check(epsilon,ifuncreal,
c    >                   bcrbot,bcibot,bcrtop,bcitop,
c    >                   ifbcbot,ifbctop,rhsinr,rhsini,
c    >                   fr,fi)
c      write(*,*) '--------------------------------'
C-Calculate analytical solution and modify existing one
C-if BC's are non-homogeneous
       do k=1,nz
         zz=z(k)
         ar(k) = analyticf(zz,zlen,alpha0)
         if (codesol(k).ne.0) then
           error(k) = (fr(k)-codesol(k))/codesol(k)
         else
           error(k) = 0.
         endif

       enddo

C-Dump out numerical and analytical solutions
       do k=1,nz
!        write(*,*) z(k),rhsinr(k),fr(k)!,fi(k)
         write(700,'(1x,5(2x,E15.8))') z(k),fr(k),ar(k),
     >   error(k),abs(error(k))
       enddo

  

       close(700)

       return
       end


C************************************************************
       function analyticf(zz,ybar,alpha0)

       pi = 4.0*atan(1.)
       a = 1./(2.*(exp(2.*ybar)+exp(-2.*ybar)))
c      a = exp(-ybar)/(exp(ybar)+exp(-ybar))
c      b = exp(ybar)/(exp(ybar)+exp(-ybar))
c      analyticf = a*exp(ybar*zz) + b*exp(-ybar*zz)
c      analyticf = sin(zz)/cos(1.)-zz
c      analyticf = sinh(2.*zz)/sinh(2.*ybar)
       analyticf = 1.
c      analyticf = 0.5*zz*zz
c      analyticf = sin(2.*pi*alpha*zz)
c      analyticf = sin(256.*pi*zz)
c      analyticf = a*(exp(2.*zz)-exp(-2.*zz))
c      analyticf = cos(2.0*pi*zz+0.5*pi)
c      analyticf = 0.5*(zz-1.0)**2.
c      analyticf = sin(0.5*pi*zz)

       return
       end
      
        
   
 

