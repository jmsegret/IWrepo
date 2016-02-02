C***************************************************
C-This file contains all routines necessary for
C-computing first two modes of solitary wave forcing
C-(PD: 3/30/04)
C***************************************************

C***************************************************
      subroutine eigenv
C***************************************************
C-Central Subroutine
C-Sets up matrices of generalized eigenvalue problem:
C-A*x = lambda*B*x
C-resulting from the Taylor-Goldstein eqn.
C-
C-x is phi. For A & B see PD's notes.
C-
C-lambad is c**2, where c is the phase speed of the
C-specific mode.
C-
C-NOTE: This subroutine is used only for a single domain
C-calculation. Data is then projected onto desired
C-multidomain grid.
C-
C-NOTE-2: nzaux = number of grid-points on multidomain grid.
C-        nzaux2 = number of grid points on uniform grid
C-                 use for intermediate calculations


      include 'dim.h'

      include 'comcoeff.h'
      include 'comleng.h'
      include 'comlegendre.h'
      include 'comsubd.h' 
      include 'comparamphys.h'
      include 'comgridptxyz.h'
      include 'comsoliton.h'

      parameter (np=nz-2,nzaux2=250)
      real amach
      real a(np,np),b(np,np) 
      complex alpha(np),evec(np,np)
      real beta(np)
C-These are the rescaled 1-2 modes desired
C-for the soliton calculation
      real test(nz),dtest(nz)
      real eigvfout(nzaux),deigvfout(nzaux),
     >d2eigvfout(nzaux)
      dimension aleg(0:nzm)
      dimension zunif(nzaux2),funif(nzaux2),
     >dfunif(nzaux2),d2funif(nzaux2) 

c     external amach,umach,dgvlrg,dgvcrg,dgpirg
 

      open(207,file='eigenv.dat')
      open(208,file='eigenf.dat')
      open(209,file='eigenf_unif.dat')
      write(207,*) 'VARIABLES = "k","Numerical"'
      write(208,*) 'VARIABLES = "k","z","`f","d`f/dz","d^2`f/dz^2"'
      write(209,*) 'VARIABLES = "k","z","`f","d`f/dz","d^2`f/dz^2"'
    
      open(210,file='test.dat')
      open(211,file='test2.dat')
      write(*,*) 'Eigenfunction Calculation'
      pi = 4.*atan(1.0)
C-Set initial guess for eigenvalue here
      ceig = brunt*zlen/(pi)
c     write(*,*) 'C=',ceig
      ushear = 0.

C-Consider a single domain and set up matrix E
      ybar = zlen

C-Set up MATRIX-A (Stratification and background current)
C
      do k1=2,nz-1
	  zz = z(k1)
          c1 = (1. - fcshear*fshear(zz,zlen)*float(isignc))**2.
        a(k1-1,k1-1) = bruntz(zz,zlen,brunt)**2./c1
      enddo 
C-Set up Matrix-B
C-(2nd order diffusive operator)

      do k1=2,nz-1
        do k2=2,nz-1

          b(k1-1,k2-1) = -((2./ybar)**2.)*d2(k1-1,k2-1)

        enddo 
      enddo


C-Call IMSL routine that returns eigenvalues (and eigenvectors)
C-(MUST CALL DOUBLE-PRECISION VERSION. OTHERWISE I GET
C-BUS ERROR MESSAGE-PD: 3/30/04).
C-
c-Additional caution: IMSL normalizes the eigenvector to
c-have Euclidean = 1. Take this into account.
c-To calculate strictly eigenvalues use routine below
c     call dgvlrg (np, a, np, b, np, alpha, beta)
c-To calculate both eigenvectors and eigenvalues
c-use following:
c-PD-2/21/05: When running on Linux machine I have
c-opted to comment out the two lines below because
c-IMSL is not installed on this machine.
c     call dgvcrg (np, a, np, b, np, alpha, beta,
c    >evec, np)

C-Calculate the performance index
c     pindx = dgpirg (np, np, a, np, b, np, alpha, beta,
c    >evec, np)
      write(*,*) 'Performance Index ',pindx

      write(*,*) 'EIGENVALUES dumped out'
      do k=1,np
        write(207,*) k,sqrt(real(alpha(k))/beta(k)) ! ,ceig/float(k)
      enddo

      cphase = sqrt(real(alpha(1))/beta(1))
      write(*,*) 'Phase Speed C=',cphase
      frn = fcshear*cphase/(brunt*zlen)
      write(*,*) 'Froude number F=',frn
      write(*,*) 'Shear amplitude=',frn*(brunt*zlen) 
C-Rescale eigenvector and transfer to (1,...,nz) array
      vmax = 0.
      do k=1,np
        var = real(evec(k,nmode))
        if (abs(var).gt.vmax) vmax = var
      enddo

      write(*,*) 'Maximum Value of Eigenvector',vmax
      eigvf(1)=0.
      eigvf(nz)=0.
      do k=1,np
        var = real(evec(k,nmode))
        eigvf(k+1) = var/vmax
      enddo
C-Interpolate eigenvector on uniform grid (z=0 to zlen)
C-Use rational interpolation routine from Numerical Recipes
      dz0 = zlen/float(nzaux2-1)
      do k=1,nzaux2
        zunif(k) = dz0*float(k-1)
        zin = zunif(k)
        call ratint(z,eigvf,nz,zin,outunif,error)
        funif(k) = outunif
      enddo

C-Compute 1st & 2nd derivative of eigenfunction on uniform grid
C-Use 2nd order finite difference 
      do k=2,nzaux2-1
        dfunif(k) = 0.5*(funif(k+1) - funif(k-1))/dz0
        d2funif(k) = (funif(k+1) - 2.*funif(k) + funif(k-1))/(dz0*dz0)
      enddo
 
      dfunif(1) = 0.5*(-funif(3) + 4.*funif(2) -3.*funif(1))/dz0
      d2funif(1) = (2.*funif(1) -5.*funif(2) +
     >4.*funif(3) -funif(4))/(dz0*dz0)

      dfunif(nzaux2) = 0.5*(funif(nzaux2-2) - 
     >4.*funif(nzaux2-1) +3.*funif(nzaux2))/dz0
      d2funif(nzaux2) = (-2.*funif(nzaux2) +5.*funif(nzaux2-1) -
     >4.*funif(nzaux2-2) +funif(nzaux2-3))/(dz0*dz0)

c     do k=1,nzaux
c       write(210,'(4(2x,E12.6))') zunif(k),funif(k),
c    >dfunif(k),d2funif(k)
c     enddo

c     close(210)

C-Calculate 1st & 2nd derivatives of eigenfunction
C-REMINDER: The calculation below is accurate only
C-when eigenfunctions are calculated for a single domain
C-calculation
C-May have to drop the following due to error in spectral
C-differentiation at boundaries !
c     ybar = zlen
c     call dpdzcol(eigvf,deigvf,d,zpts,nzm,ybar,1,nzloc)
      
c     call d2pdz2col(eigvf,d2eigvf,d2,zpts,nzm,ybar,1,nzloc)

c     do k=1,nz
c       write(211,'(4(2x,E12.6))') z(k),eigvf(k),
c    >deigvf(k),d2eigvf(k)
c     enddo      
C-output eigenvector
      write(*,*) 'EIGENVECTOR dumped out'
      sum = 0. 
      sumsol = 0.

C-Final interpolation of eigenfunction from uniform grid
C-to non-uniform multidomain grid
      do k=1,nzaux
        zin = zaux(k)
        call ratint(zunif,funif,nzaux2,zin,outunif,error)
        eigvfout(k) = outunif
        call ratint(zunif,dfunif,nzaux2,zin,outunif,error)
        deigvfout(k) = outunif
        call ratint(zunif,d2funif,nzaux2,zin,outunif,error)
        d2eigvfout(k) = outunif
      enddo

      do k=1,nzaux
        write(208,'(1x,I3,4(2x,E12.6))') k,zaux(k),
     >eigvfout(k),deigvfout(k),d2eigvfout(k)

        write(209,'(1x,I3,4(2x,E12.6))') k,zunif(k),
     >funif(k),dfunif(k),d2funif(k)

      enddo

      close(207)
      close(208)
 
      
      return
      end

C

C***************************************************
      subroutine kdv_coeff(alphak)
C***************************************************
C-This subroutine calculates the coefficients
C-of the non-dimensional KdV equation.
C-Input is simply eigenfunction.
      include 'dim.h'

      include 'comcoeff.h'
      include 'comleng.h'
      include 'comsubd.h'
      include 'comparamphys.h'
      include 'comgridptxyz.h'
      include 'comsoliton.h'

      real intk

c     dimension eigvf(nz)
      dimension fint(nz) 
      dimension f(nzloc),df(nzloc),a(nzloc)

      open(308,file='test_unix.dat')

C-Note: All coefficients could be calculated
C-simultaneously but I prefer to play it conservative
C-(PD: 4/5/04).
C
C-Calculate Integral I_k

      intk = 0.

      do ks=1,nsubd
        ybar = zh(ks)
      
        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          zz = z(k)
          fnum = bruntz(zz,zlen,brunt)**2. * 
     >           fsol(k)**2.
          fden = (umax**3.)*(1. - 
     >           fcshear*fshear(zz,zlen)*float(isignc) )**3.
          a(kloc) = fnum/fden 
c         write(*,*) zz,a(kloc) ! ,fnum,fden
        enddo

        call integrcol(a,outint,nzm,ybar)
c       write(*,*) nzm,ybar,outint
        intk = intk + outint

      enddo

      cnorm = -2.*brunt*zlen*zlen
      intk = intk * cnorm

      write(*,*) 'I_k ',intk

C-Calculate beta_k
      betak = 0.

      do ks=1,nsubd
        ybar = zh(ks)

        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          a(kloc) = fsol(k)**2.
        enddo

        call integrcol(a,outint,ndim,ybar)

        betak = betak + outint

      enddo

      cnorm = -1./(intk*zlen) * (brunt*zlen**3.)
      betak = betak * cnorm

      write(*,*) 'beta_k',betak ! /(zlen**3.*brunt)

C-Calculate alpha_k
      alphak = 0.

      do ks=1,nsubd
        ybar = zh(ks)

C-Need to first compute derivative of N^2(z)
        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          zz = z(k)
          f(kloc) = bruntz(zz,zlen,brunt)**2. 
        enddo

        call dpdzcol(f,df,d,zpts,nzm,ybar,1,nzloc)

        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          zz=z(k)

          fnum1 = 3.*fcshear*f(kloc)*dfsheardz(zz,zlen)*float(isignc)
          fden1= (1.-fcshear*fshear(zz,zlen)*float(isignc))**4.

          fnum2 = 2.*df(kloc)
          fden2 = (1.-fcshear*fshear(zz,zlen)*float(isignc))**3.
         
c         write(308,*) zz,3.*fcshear*f(kloc)*1.*float(isignc)/
c    >(1.-fcshear*fshear(zz,zlen)*float(isignc)),
c    > (2.*df(kloc)) !fnum1,fden1
c         write(*,*) zz,f(kloc),df(kloc) 
          a(kloc) = (fnum1/fden1+fnum2/fden2)*(fsol(k)/umax)**3.
c         a(kloc) = fnum2*(fsol(k)/umax)**3.
        enddo

        call integrcol(a,outint,ndim,ybar)

        alphak = alphak + outint

      enddo

      cnorm = brunt*(zlen**3.)/(intk*zlen)
      alphak = alphak * cnorm ! / zlen

      write(*,*) 'alpha_k',alphak

C-Now calculate Soliton wavenumber (k_s)
C-CAREFUL: When running shear case, the definition of alpha0
C-should be changed !
      alpha0 = (zeta0/zlen)*abs(umax)*zlen
c     alpha0 = (zeta0/zlen)*abs(brunt*zlen)*zlen

      solks =  sqrt(alpha0*abs(alphak/(betak*12.)))
c     solks = 0.7180
      write(*,*) 'Soliton Wavenumber',alpha0,solks,1./solks
        
      return
      end
      

C****************************************************
      function strat(z,ceig,brunt,ushear)
C****************************************************
C-Coefficient of stratification term in Taylor-Goldstein
C-Equation.
C-PD( 3/30/04): IMPROVE THIS !

      freq2 = bruntz(z,zlen,brunt)**2.
      strat = freq2/(ushear - ceig)

      return
      end

C****************************************************
      function tgsol(z,zlen)
C****************************************************
C-Analytical solution of TG eqn. for uniform stratifcation
C-(Test case only)
      pi = 4.*atan(1.0)

      tgsol = sin(2.*pi*z/zlen)
c     tgsol = exp(0.5*z*z)

      return
      end
