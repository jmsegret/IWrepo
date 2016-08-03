c^^^^^^^^^^^^^^^^^^^^^^^^^^^
      Program Dubwise 
c^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c  A Three-Dimensional Unsteady Navier-Stokes Solver
c
c  copyright (c) Peter Diamessis, 2002 
c  The original code was developed in 1993 by Dan Chan.
c  This code employs a spectral subdomain decomposition method
c  in the vertical combined with spectral filtering when higher
c  Reynolds numbers are used
c
c  homogeneous in x and y direction
c  inhomogeneous in z direction
c
c  specified temperature at the boundaries
c  no slip velocity at the lower boundary (k=1)
c  stress free upper boundary (k=nz)
c  The boundary conditions can easily be modified in the new
c  code structure.
c  
c  //////////////////////////////////////////////
c  output send to t_(time) every iform step
c  by plotf (u,v,w,t) in physical space
c
c  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C-PD: 2/24/04: Hereafter the variable temp actually
C represents the density (Not temperature) perturbation
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



c///////////////////////////////////////////////////
      character*20 fout
      character*7 fmt1

c//////////////////////////////////////////////////
      include 'dim.h'

c 
c//////////////////////////////////////////
      include 'paramvarfilterleg.h'
      parameter (imaxcalc=1)
c-iforcelfow: Run a forced flow simulation without any check for relaxation.
c-irelax: Perform relaxation process.
c-irelax: Flag indicating whether relaxation process has reached end.
c /////////////////////////////////////////
c
      include 'comwave.h' 
      include 'comleng.h'
      include 'comtime.h'
      include 'comflow.h'
      include 'comfft.h'
      include 'commatrix.h'
      include 'comcoeff.h'
      include 'comparamphys.h'
      include 'combc.h'
      include 'comscratch.h'
      include 'commass.h'
      include 'comnusselt.h'
      include 'commeanv.h'
      include 'comzeta.h'
      include 'combcg.h'
      include 'comsubd.h'
      include 'comspecfilter.h'
      include 'compencoeff.h'
      include 'comtimeadvance.h'
      include 'comsoliton.h'
      include 'comgridptxyz.h'
      include 'comcritlayer.h'

c************************************

c     dimension apvel(nz,nz,nxhpy),app(nz,nz,nxhpy),apt(nz,nz,nxhpy),
c    1          apw(nz,nz,nxhpy)
c     dimension kpv(nz,nxhpy),kpp(nz,nxhpy),kpt(nz,nxhpy),
c    1          kpw(nz,nxhpy)

      dimension xdudz(nxpp,ny,nz),xdwdx(nxpp,ny,nz),xdtdz(nxpp,ny,nz)


      dimension to(17)

      character*2 f11
      character*12 fil
      character*9 f12
      character*8 f13
      character*25, fil1


      real nuss

      logical lwrite,cpenspecfilter,ldtadapt

c-IMSL related bullshit for timing
      integer isec0,isec1
      integer isecs,ihourf,iminf,isecf
      real timediff,timearray(2)



       common /deriv/ xdudz,xdwdx,xdtdz

       common /Dump/ Tdump,Tx     
       common /anim/ianim,ianimc
       common /Misc/tzmode
c       common  /forcef/  phiz(nz)
       common /domlength/ ZL,Zcur

c          tzmode=100000.0        
 
           fcshear=0.000

c
c-------------------------------------------------------------------
c  KEY PARAMETERS (INPUT etc.)
c-------------------------------------------------------------------
c
c  istart:  set to 1 to restart from a previous run
c  nsteps:  number of time steps to run
c  dt:      size of time step
c  nwrt:    no. of timesteps one dumps out sample z-profile
c  iform:   format of 3-D field output (0 = binary)
c
c  xlen:    domain length in x-direction.
c  ylen:    domain length in y-direction
c  zlen:    domain length in z-direction
c  alpha:   wavelength in the x-direction
c  beta:    wavelength in the y-direction
c  brunt:   Reference value for Brunt Vaisala Frequency
C           (now an external input)
c  pra:     Actual  Prandtl number (no longer inputted)
c  raa:     Actual Rayleigh  number
c  ta:      Taylor number (not used)
c  xnu:     kinematic viscosity
c  xkappa:  thermal diffusivity
c  grav:    gravitational acceleration
c  expan:   thermal expansion coefficient
c  omega:   ambient rotation rate.
c  umax:    maximum reference velocity 
c           (centerline velocity for wake/
c            mean stream velocity for resuspension problem)
c  tempgrad: Non-dimensional value of ambient temperature gradient.
c  rho0: reference density 
c  zeta0:    soliton amplitude (resuspension only)
c  fcshear:  U_0/c (Ratio of Max. Shear Velocity to phase speed)
c  isignc: sign of wave (upstream propagation = -1)
c  nmode:    Wave mode #.
c
c  iforceflow: forcing turned off/on (no longer used)
c  irelax: is relaxation used ?
c  ispecfilternl: use filtering after N-L term calculation
c  ispecfilterexp: use filtering during pressure and viscous term
c                  calculation
c  nspecfilterexp: frequency of pressure/viscous solution filtering
c  amp:       exponential filter amplitude
c  p:         2 x order of exp. filter (Legendre space).
c  pf:        2 x order of exp. filter (Fourier space)
c  fackc: Lag of exp. filter
c  fackcbv: Lag of Boyd-Vandeven filter
c  facdir: Factor multiplying coefficient in Dirichlet BC
c  facneummn: -''-        -""-         -""- in Neumann BC
c  facrobin:  -''-        -""-         -""- in Robin BC (patching condition)
c  navgintrfc: No. of timesteps to perform interfacial averaging.
c  ncalclegspec: No. of timesteps to calculate Legendre spectra.
c  ifilter: Apply estimation filter (y/n) ?
c  nstepfilter: No. of timesteps to apply estimation filter.
c  tstart: Physical time (sec) at which is started (this is used
c          only when restarting with initial data from relaxation)
c  itempreset: Reset temperature perturbations upon restart (y/n)
c  nrestartdump: No. of timesteps to dump out restart file.
c  n3doutput: No. of timesteps to dump out 3-D (u,v,w,rho) field.
c  i0 = 0: For given value of dt, when restarting a run, how many
c          timesteps from t=0 is restart time.
c
c-------------------------------------------------------------------
      namelist/inputs/ istart,dtin,nsteps,tanim,iform,cdtin,waveperiods
      namelist/params/ xlen,ylen,zlen,brunt,xnu,xkappa,grav,
     >expan,umax,rho0,xacrit,xbcrit,zcen,xacu,zsh

      namelist/mflags/irelax,ispecfilternl,
     >ispecfilterexp,nspecfilterexp,amp,p,pf,fackc,fackcbv,
     >facdir,facneumn,facrobin,navgintrfc,ncalclegspec,
     >ifilter,nstepfilter,
     >tstart,npostp,itempreset,nrestartdump,n3doutput,i0

          
c       real pi2,tshutoff,xacrit,xkcrit,xmcrit,xomega,uo,zcen,xbcrit,zdim

     
c      
C-PD-2/20/02: COMMENT THIS LINE OUT AND ENTER YBAR THROUGH INPUT FILE
C-YBAR: Height of domain.
c     ybar=1.0
c
c  set up the working arrays for the FFT calls
c
c     call horfft_init
      call fftfax(nx,ifaxx,trigsx)
!-If this is a 3-D run, set up working arrays for y-FFT calls
      if (ndims.eq.3) then
        call cftfax(ny,ifaxy,trigsy)
!-Call below is used only when working with variables
!-that are funcitons of (y,z) (e.g. mean fields).
!-Allows real-to-complex FFT's in y-direction.
        call fftfax(ny,ifaxy1d,trigsy1d) 
      endif
c     call cftfax(nz,ifaxz,trigsz)

      scalez= 1./float(nz)
      scaley = 1./float(ny)
c

      pi2 = 8.*atan(1.0)

c
      do i=1,nxpp
        do j=1,ny
          do k=1,nz
           
          su(i,j,k) = 0.
          sv(i,j,k) = 0.
          sw(i,j,k) = 0. 
          st(i,j,k) = 0.
            
          enddo
        enddo
       enddo
c
c-The following has been commented out as BC's for velocity are set
c-in init
c--- specify boundary conditions for velocity at initial time - will
c     be modified if restart file is used later
c     first 1,nxppy elements for the lower wall
c     nxppy+1,2*nxppy for the upper wall
      
c     do i=1,nxpp
c       do j=1,ny
c         do k=1,2
c          uw(i,j,k) = 0.
c          vw(i,j,k) = 0.
c          ww(i,j,k) = 0.
c         enddo
c       enddo
c     enddo
c---  end specification of b.c. for velocity
c     boundary conditions for temperature in subroutine init.f
c     
 
      open(7,file='input.dat',status='old')
c     open(10,file='time.dat',status='unknown')
      open(9,file='nu.dat',status='unknown')
c     
      read(7,inputs)
      read(7,params)
      read(7,mflags)
c
C-Output all parameters as reminder
      write(*,*) '*-------------------------------*'
      write(6,inputs)
      write(6,params)
      write(6,mflags)
      write(*,*) '*-------------------------------*'

C-OPEN OUTPUT FILES HERE
      open(50,file='minmaxtest.dat',status='unknown')
      open(51,file='cfl.dat',status='unknown')
      open(55,file='tkepe.dat',status='unknown')
      open(60,file='energy.dat',status='unknown')
      open(61,file='energy_wake.dat',status='unknown')
      open(62,file='wakescale.dat',status='unknown')
      open(63,file='avgint.dat',status='unknown')
      open(64,file='limitcycle.dat',status='unknown')
      open(65,file='dtchange.dat')
      open(66,file='growthrate.dat')
      open(67,file='animdump.log')
 
C-OUTPUT FILE HEADERS
      write(51,*) 'VARIABLES = "CFLXMAX","I","J","K",
     >                         "CFLYMAX","I","J","K",
     >                         "CFLZMAX","I","J","K"'
      write( 55,*) 'VARIABLES = 
     >"TIME","Nt","TKE","PRODY","PRODZ","PROD","BFLUX",
     >"EPSILON","EPSSGS","EPSTOT",
     >"P/`e_t"'
c
      write(60,*) 'VARIABLES =
     >"TIME","Nt","TKE","u1tke","u2tke","u3tke","PRODY","PRODZ",
     >"PROD","BFLUX","CHI","CHISGS","CHITOT",
     >"EPSILON","EPSSGS","EPSTOT",
     >"TKEmean","Tvar"'
c
      write(61,*) 'VARIABLES =
     >"TIME","Nt","TKE","u1tke","u2tke","u3tke","PRODY","PRODZ",
     >"PROD","BFLUX","CHI","CHISGS","CHITOT",
     >"EPSILON","EPSSGS","EPSTOT",
     >"TKEmean","Tvar"'
      write(62,*) 'VARIABLES =
     >"TIME","Nt",
     >"umax","umaxy","umaxz","L_y","L_z","yc","zc"'
      write(63,*) 'VARIABLES = "time","pctu","pctv","pctw","pctt"'
 
      write(64,*) 'VARIABLES = "t","u_1","w_1","u_2","w_2","u_3","w_3"'

      write(66,*) 'VARIABLES = "t","wrms","wmax","zmax"'

      close(66)

      close(67)
c    ZL is used in phiz inside forcing.f 
   
          ZL=zlen      
          Zcur=zsh


c ----begin operating conditions and material properties---                                                                      
      alpha = pi2/xlen
      beta  = pi2/ylen
      gamma = pi2/zlen

!----------------------------------------------------------
!-Define wave packet and localization function values here
!----------------------------------------------------------
!AMA 02-28-2008 to improve vertical resolution Lambdax/Lambdaz was limited to 4
! Lambdax=xlen and the length scale was chosen as Lambda x similar to Winter's and D'Asaro
!-PD-Lyon 7/1/15: Modified Ammar's approach to work with MODE-1 wave
!
       zdim= xlen ! -PD-7/1/15: Horizontal wavelength
                  !             following TCFD12 paper
       xkcrit=pi2/xlen !-Horizontal wavenumber
       xmcrit=pi2/(2.*zlen) !-Vertical wavenumber
                            ! CAREFUL: Domain is designed to
                            ! accommodate half a wavelength
       xomega=(bruntz(0.0,zlen,brunt)*xkcrit)
     >           /sqrt(xkcrit**2+xmcrit**2)      !-Wave Frequency
       uo= 0.07

       cphx = xomega/xkcrit !-Wave horizontal phase speed

!-Output Background &  Wave Characteristics
!
      write(6,*) '---------------------------------------------'
      write(6,*) 'WAVE AND BACKGROUND CHARACTERISTICS:'
      write(6,*) 
      write(6,*) 'B.V. Frequency (rad/s): ', bruntz(0.0,zlen,brunt)
      write(6,*) 'B.V. Period (s): ', pi2/brunt
      write(6,*) 'Mean density gradient: ',
     >            rhograd(z,zlen,grav,rho0,brunt),
     >            brunt*brunt*rho0/grav
      write(6,*) 'Domain length and height (m): ',xlen,zlen
      write(6,*) 'Number of grid points in x and z dirs.: ',
     >            nx,nzloc*nsubd
      write(6,*) 
      write(6,*) 'Wave Frequency omega (rad/s): ',xomega
      write(6,*) 'Wave Period (s): ', pi2/xomega
      write(6,*) 'Characteristic time lambda_x/U0 (s): ', xlen/umax
      write(6,*) 'Max. horiz. wave current U0 (m/s): ',umax
      write(6,*) 'Forcing amplitude function: ', xacrit
      write(6,*) 'Horizontal Phase Speed (m/s): ',cphx
      write(6,*) 'Horizontal and vertical wavenumbers (rad/m)',
     >           xkcrit,xmcrit
      write(6,*) '---------------------------------------------'

       open(543,file='waveparams')
       write(543,*)'omega=',xomega,'rad/s T=',pi2/xomega,'sec Cph=',
     > xomega/sqrt(xkcrit**2+xmcrit**2),'m/s Cg=',
     > xomega*xmcrit/(xkcrit*sqrt(xkcrit**2+xmcrit**2)),'m/s','Amp='
     > ,xacrit,uo
c        open(5786,file='testvars')
c        write(5786,*)'pi2=',pi2,'xacrit,xkcrit,xmcrit,xomega,uo,zcen,
c     >  xbcrit,zdim,xlen',xacrit,xkcrit,xmcrit,xomega,uo,zcen,
c     >  xbcrit,zdim,xlen
      if(istart.eq.0) then

      open(698,file='xtu.dat')
      write(698,*)'VARIABLES ="x","time","u","u_mid","u_pyc"'
      write(698,*)'ZONE I=',nx,',J=',1+(nsteps/4),',',' F=POINT'
      close(698)
      open(697,file='xtw.dat')
      write(697,*)'VARIABLES ="x","time","w","w_mid","w_pyc"'
      write(697,*)'ZONE I=',nx,',J=',1+(nsteps/4),',',' F=POINT'
      close(697)
      open(696,file='xttemp.dat')
      write(696,*)'VARIABLES="x","time","delN2","delN2_mid","delN2_pyc"'
      write(696,*)'ZONE I=',nx,',J=',1+(nsteps/4),',',' F=POINT'
      close(696)

      open(700,file='uveltest.dat')
      close(700)
      else
      continue
      endif

c nondimensional rotation rate omega*h^2/xkappa      
c-(**PD_modify**: Introduce externally specified 
c-dimensional rotation rate (Coriolis
c-parameter ?), which is not used right now)
c      omz = 0.5*pr*sqrt(ta)
       omz = omega

c
c thermal conductivity
c(**PD_modify**: Thermal diffusivity now dimensional.
c-Read in from input file)
c      xkappa=xnu/pr
c-Set the actual Prandtl #:
      pra = xnu/xkappa

c
c-Set Prandtl number in code (In reality this is the kin. viscosity nu)
      pr=xnu
C-Set Relaxation constant to zero
      irelaxend = 0
c-Set rayleigh number in code (This is not the actual Rayleigh #).
c-It is assumed that domain height is equal to H=1 here.
c      ra=xkappa*raa

c initial temperature drop across the layer (dimensional in K !)
c used only as a temperature scale to nondimensionalize equations
c
C(**PD_modify**: Now use actual Rayleigh #)
c     tdiff = xkappa*xnu*raa/(expan*grav)
C-PD: 2/24/04. Mean Density Gradient calculated in setup now

c  buoyancy frequency based on tdiff
c     brunt=sqrt(xkappa*xnu*raa/(zlen**4.))                  
C-PD: 2/24/04. BV-Freq. supplied externally
c

c  save initial spectral filter order
c  then set to value needed for first step of gradual filter increase
      pinit = p
c     p = 3.
c-Interval for increasing wave amplitude
      nbump = 50
 
      print *,'xlen, ylen,zlen:     ',xlen,ylen,zlen
      print *,'alpha, beta:    ',alpha,beta
c     print *,'Ra,Ta,Pr:       ',ra,ta,pr
      print *,'Ra,Ta,Pr:       ',raa,ta,pra 
      print *,'nu,kappa,expan: ',xnu,xkappa,expan
      print *,'delT,freq,omega:',tdiff,brunt,omz
      print *,'iforceflow,irelax,ifilter,nstepfilter,tstart,npostp,
     +itempreset,nrestartdump,n3doutput',
     +iforceflow,irelax,ifilter,nstepfilter,tstart,npostp,itempreset,
     +nrestartdump,n3doutput

C-Set everything relevant to timestep here
      dt = dtin
      cdt = cdtin

C-Set factors in Helmholtz solvers right here (They are now needed
C-for setup)
      fact = 1./(xkappa*dt)
      facvel = 1./(pr*dt)

c-Set Output Times
c-PD: Refine this later
        to(1) = 400.
	to(2) = 970.
	to(3) = 1370.
	to(4) = 1770.
	to(5) = 2170.
	to(6) = 2570.
	to(7) = 3174.
	to(8) = 3655.
      do i=9,17
        to(i) = 4000. + float(i-9)*200.
      enddo

c-Set animation flag & counter to zero
      ianim = 0  
c      ianimc = 0

c      
c ---- end operating conditions and material properties---  
c     
      if(istart.eq.0) then
        call timecall(isec0)
        write(*,*) 'SETUP'
        call setup
        call timecall(isec1)
        call timeelapse(isec0,isec1)

        write(*,*) 'MALAKAS', xnu,xkappa
 
!       write(*,*) 'Calculating Forcing Eigenfunctions'
!       call eigenv 
!       stop
               

        write(*,*) 'INITIALIZATION'
        call timecall(isec0)
        call init
!       call testdump(t,ybar,0,tstart)
c       isgscalc = 0
c       call filtering(isgscalc)
        call timecall(isec1)
        call timeelapse(isec0,isec1)

        t = 0.

C-Update timestep file      
c        write(65,'(A10,F15.8,A21,F8.3)') 'At time t=',t,
c     > ', timestep set to Dt=',dt
c-For INITIAL CONDITION
c-Calculate Legendre spectra of u,v,w,T at y=ny/2
c      write(*,*) 'CALCULATING LEGENDRE SPECTRA of I.C.'
c      call speclegencalc('u',t,u)
c      call speclegencalc('v',t,v)
c      call speclegencalc('w',t,w)
c      call speclegencalc('T',t,temp)

c-output initial u,v,w,T field
        f11 = 'p_'
        write(fil,15)f11,t-tstart

        if (iform.eq.0) then
          open(3,file=fil,form='unformatted',status='unknown')
        endif

        if (iform.eq.1) then
          open(3,file=fil,status='unknown')
        endif

        call plotf(3,iform,t)
        close(3)

        fac1 = 1.0
        fac2 = 0.
        fac3 = 0.
        bdf1 = 1.
        bdf2 = 0.
        bdf3 = 0.
        bdfv = 1. 
        initenergy = 0

      else 

        open(8,file='start.dat',form='unformatted',status='unknown')
        call restart(8,t) 
        close(8)

C-CAREFUL: This is temporary because ra needs to be reset
C-Startup value is 0.
c        ra=xkappa*raa

c-Reset initial temperature field
c-This is done when we need to start a run with the
c-end result of a relaxation rpocess. T must only have mean
c-component at t=0.
        if (itempreset.eq.1) call tempreset(zlen)
!-PD-10/13/08: The following line has been uncommented to
!-recompute wave#'s, grid-points, diff. matrices when restarting
!-on a new grid

      call setup
      write(*,*) 'MALAKAS', xnu,xkappa


c-output initial u,v,w,T field
        f11 = 'p_'
        write(fil,15)f11,t-tstart

        if (iform.eq.0) then
          open(3,file=fil,form='unformatted',status='unknown')
        endif

        if (iform.eq.1) then
          open(3,file=fil,status='unknown')
        endif

        call plotf(3,iform,t)
        close(3)
C-Update timestep file
        write(65,'(A10,F15.8,A21,F8.4)') 'At time t=',t,
     > ', timestep set to Dt=',dt

C-PD: 10/9/03: This is a very special modification that
C-will allow restart with variable timestep to maintain
C-3rd order accuracy in BDF approximation of non-linear
C-term. The following coefficients have to do with the very
C-first step after restarting. We're interested in a doubled
C-timestep but this may be modified to be more general
        dt1 = dt
        dt2 = cdt*dt
        dt3 = cdt*dt

C-AB3-Variable timestep
        fac1 = abvar1(dt1,dt2,dt3)
        fac2 = abvar2(dt1,dt2,dt3)
        fac3 = abvar3(dt1,dt2,dt3)

C-BDF3-Variable timestep
        bdf1 = bdfvar1(dt1,dt2,dt3)
        bdf2 = bdfvar2(dt1,dt2,dt3)
        bdf3 = bdfvar3(dt1,dt2,dt3)
        bdfv = bdfvarv(dt1,dt2,dt3)

        initenergy = 1
        call postp(initenergy,t,brunt,tstart)

!       call testdump(t,ybar,0,tstart)
      end if

      lnoise = .false.

c     call testdump(t,ybar,0,tstart)

C-Begin main solution loop here
      call timecall(isecs)
c     timediff = dtime(timearray)
      do 1000 istep = 1,nsteps

c-Calculate minimum and maximum values of (u,v,w,T) and identify
c-their k-index position
        if (imaxcalc.eq.1) call maxcalc(istep)

c-If adaptive timestepping is turned on advance IDTSTEP flag
        if (ivardt.eq.1) then
          idtstep = idtstep+1
c          write(65,*) 'Global timestep, Adaptive timestep',
c     >    istep,idtstep
        endif

        write(*,*) '*************************************'
        write(*,*) 'Now Solving Step',istep
        write(*,*) '*************************************'
c




c           if((t.gt.5).and.lnoise) then
c           call addVelocityNoise(istep,xlen)
c           lnoise = .false.
c           end if
c


c-Calculate convective terms. RHSFM remains unchanged from
c-2nd order NS solver.
        write(*,*) 'CALCULATING CONVECTIVE TERMS'
        call timecall(isec0)
        call convec(istep,istart,ispecfilternl)
        call timecall(isec1)
        call timeelapse(isec0,isec1)





c
c-Adjust Non-Linear term advancement coefficients 
c-(either stiffly stable or AB) and temporal derivative
c-(everything except with the one associated viscous
c-equation).

C-Am I still timestep at n<2 (start or restart regardless) ?
        if (istep.lt.2) then

C-Start from scratch, use SS2/BDF2
          if(istart.eq.0) then
C-AB2
c           fac1 = 1.5
c           fac2 = -0.5
c           fac3 = 0.
C-SS2
            fac1 = 2.
            fac2 = -1.
            fac3 = 0.
C-BDF2
            bdf1 = 2.
            bdf2 = -0.5
            bdf3 = 0.
C-Otherwise restarting. Use AB3/BDF3 with variable timestep
C-flexibility

          else

            dt1 = dt
            dt2 = dt
            dt3 = cdt*dt

C-AB3-Variable timestep
            fac1 = abvar1(dt1,dt2,dt3)
            fac2 = abvar2(dt1,dt2,dt3)
            fac3 = abvar3(dt1,dt2,dt3)

C-BDF3-Variable timestep
            bdf1 = bdfvar1(dt1,dt2,dt3)
            bdf2 = bdfvar2(dt1,dt2,dt3)
            bdf3 = bdfvar3(dt1,dt2,dt3)

          endif

        else

C-Am I at timestep n=2 or beyond ?

          ldtadapt = ((ivardt.eq.1).and.(idtstep.lt.2))
C-If I'm not using adaptive timestepping (start or restart
C-regardless) use BDF3/SS3.

          if (.not.ldtadapt) then
C-AB3
c         fac1 = 23./12.
c         fac2 = -8./6.
c         fac3 = 5./12.
C-SS3
            fac1 = 3.
            fac2 = -3.
            fac3 = 1.
C-BDF3
            bdf1 = 3.
            bdf2 = -1.5
            bdf3 = 1./3.

        else
C-The only other possibility is having adaptive timestepping
C-and being at timestep idtstep=1
            dt1 = dt
            dt2 = dt
            dt3 = cdt*dt

C-AB3-Variable timestep
            fac1 = abvar1(dt1,dt2,dt3)
            fac2 = abvar2(dt1,dt2,dt3)
            fac3 = abvar3(dt1,dt2,dt3)

C-BDF3-Variable timestep
            bdf1 = bdfvar1(dt1,dt2,dt3)
            bdf2 = bdfvar2(dt1,dt2,dt3)
            bdf3 = bdfvar3(dt1,dt2,dt3)
        endif

      endif
c

c
c-Below is the input for the pressure eqn. which is solved first
c-in the subsequent fractional step.
c

        do i=1,nxpp
         do j=1,ny
          do k=1,nz
            su(i,j,k)  = u(i,j,k)
            sv(i,j,k)  = v(i,j,k)
            sw(i,j,k)  = w(i,j,k)
c-In the previous code the following line was included in helmv.
            st(i,j,k)  = temp(i,j,k)
          enddo
         enddo
        enddo


C-Now calculate each primitive variable
C-First determine whether spectral filtering is needed
C-in penalty solutions.
       cpenspecfilter = ((ispecfilterexp.eq.1).and.
     > (mod(istep,nspecfilterexp).eq.0))
       write(*,*) 'SPECTRAL FILTERING of penalty solutions',
     >cpenspecfilter
c
       write(*,*) 'SOLVING VISCOUS EQUATIONS'
       call timecall(isec0)
       call calcuvwt(izmode,ispecfilter,irelax,
     > cpenspecfilter)
       call timecall(isec1)
       call timeelapse(isec0,isec1)

C-Viscous component of BDF
C-Am I at timestep < 2 (start/restart regardless) ?
        if (istep.lt.2) then

c-If starting from scratch use 2nd order BDF
          if (istart.eq.0) then
            bdfv = 1.5
          else
c-If restarting use variable timestep 3rd order BDF
            dt1 = dt
            dt2 = dt
            dt3 = cdt*dt
            bdfv = bdfvarv(dt1,dt2,dt3)
          endif

        else

C-If we are at timestep n=2 or beyond
C-(start/restart regardless)
          ldtadapt = ((ivardt.eq.1).and.(idtstep.lt.2))

C-If not using adaptive timestepping use 3rd order BDF
          if (.not.ldtadapt) then
            bdfv = 11./6.

          else

C-If using adaptive timestepping use 3rd order BDF
C-for variable timestep. Set cdt to 1 to make sure
C-no further timestep alterations are performed.
C-Given that this will be used at timestep n=2
C-of adaptive timestepping set ivardt to 0 to ensure
C-that regular timestepping will be set up henceforth.

            dt1 = dt
            dt2 = dt
            dt3 = cdt*dt
            bdfv = bdfvarv(dt1,dt2,dt3)
            cdt = 1.
            ivardt = 0

          endif

        endif

C-Gradually increase soliton wavenumber if desired.
C-Or gradually increase filter order
C-(PD: 1/24/04: Specialized version where wave# is
C-increased and updated by velocity perturbation.)
c-Call Senthil's filtering routine to attenuate higher
c-wavenumber activity
c       if ( ((zeta0.lt.0.25).and.
c     >     (zeta0.gt.0.)) .and.
c     >     (mod(istep,nbump).eq.0) )
c     >      then
c         write(*,*) 'BUMPING ZETA_0 TO',zeta0+0.025
c         zeta0 = zeta0 + 0.025
c        write(*,*) 'Resetting Wave and Boundary Conditions'
c         call wavesetup
c         call bcreset
c      endif

c       if ((mod(istep,nbump).eq.0).and.(p.lt.pinit)) then
c         write(*,*) 'INCREASING EXPONENTIAL FILTER ORDER TO',p+0.5
c         p = p + 0.5
c       endif

       if ((ifilter.eq.1).and.
     +   (mod(istep,nstepfilter).eq.0)) then
c    +   (int(to(3)/dt).eq.istep))) then
         write(*,*) 'APPLYING ESTIMATION MODEL FILTER'

C-Note the flag isgscalc is a flag which indicates whether the
C-SGS stress and heat flux tensors should be calculated immediately
C-after filtering
         if (mod(istep,npostp).eq.0) then
           isgscalc = 1
         else
           isgscalc = 0
         endif

         call timecall(isec0)
         call filtering(isgscalc)
         call timecall(isec1)
         call timeelapse(isec0,isec1)

       endif

c-Update time
       t = t + dt

c-Calculate Legendre spectra of u,v,w,T at y=ny/2
       if (mod(istep,ncalclegspec).eq.0) then
         write(*,*) 'Calculating LEGENDRE spectra'
         call speclegencalc('u',t,u)
         call speclegencalc('v',t,v)
         call speclegencalc('w',t,w)
         call speclegencalc('T',t,temp)
       endif

C-Implement Relaxation Procedure ?
      if (irelax.eq.1)
     +call relaxation(initenergy,t,brunt,
     +irelaxend)
C-Is relaxation complete ? Then simply output data to run.dat file
      if (irelaxend.eq.1) write(*,*) 'RELAXATION COMPLETED'

      if (mod(istep,navgintrfc).eq.0) then
         write(*,*) 'ADAPTIVE AVERAGING u,v,w, T
     >at subdomain interfaces'
         call interface_avg3d(u,icu)
c        pause
         call interface_avg3d(v,icv)
         call interface_avg3d(w,icw)
         call interface_avg3d(temp,ict)
         totint = float(ny*nxpp*(nsubd-1))
         write(63,'(5(1x,E12.6))') t,
     >float(icu)/totint,float(icv)/totint,
     >float(icw)/totint,float(ict)/totint
       endif

C-Is Postprocessing needed ?
      if (mod(istep,npostp).eq.0)
     >call postp(initenergy,t,brunt,tstart)

        do i=1,nxpp
          do j=1,ny
            do k=1,nz
              su(i,j,k) = 0.
            enddo
          enddo
        enddo

C-If this is an appropriate time, dump out 3-D field
C-and other parameters to disk.
C-i0 = parameter relevant for restart runs.
C-If this is a restart run with tstart<>0 and an additional
C-restart is now being performed, i0 indicates the number
C-of timesteps (of the current value of Dt) from tstart
C-until the time at which current run was restarted.
c     i0 =  1050  ! 300 ! int(tstart/dt)
      lwrite = .false.
      do m=1,17
        iistep = int((to(m))/dt)
c        lwrite = (iistep.eq.(istep+i0))
         lwrite = (t.eq.to(m))
        if (lwrite) goto 90
      enddo

C-Dumps out output file (FOR RESTART PURPOSES) at regular intervals
      if (mod(istep,nrestartdump).eq.0) then
         write (*,*) 'DUMPING OUT RESTART FILE'
c        f13 = 'run.dat.'
c        fout = f13
C-To rewrite over old restart file set file name to start.dat
c        write(unit=fout(7:9),fmt=fmt1) nint(time(it))
c        write(fout,25) f13,t-tstart
         fout = 'start.dat'
         open(2,file=fout,form='unformatted',status='unknown')
         call output(2,t)
         close(2)
      endif


c     if (istep.eq.2000000) then
c 90   if(lwrite) then
  90   if(mod(istep,n3doutput).eq.0) then
c     write(*,*) 'Dumping out 3-D data file'
c     fmt1='(i1)'
c     timeout=float(istep)*dt
c     if (nint(timeout.ge.10.) fmt1='(i2)'
c     f11 = 'run.dat.'
c     fout = f11
c     write(unit=fout(7:9),fmt=fmt1) nint(time(it))
c     open(2,file='run.dat',form='unformatted',status='unknown')
c     call output(2,t)
c     close(2)
c


      f11 = 'x_'
      f12 = '_uvwt.mat'

      write(fil1,151)f11,t,f12
 151  format(a2,e10.4,a9)


      f11 = 'p_'
      write(fil,15)f11,t-tstart
c15    format(a2,e10.4)
c
c ///////////////////////////////////
      if (iform.eq.0) then
         open(3,file=fil,form='unformatted',status='unknown')
      endif

      if (iform.eq.1) then
         open(3,file=fil,status='unknown')
      endif

      if (iform.eq.10) then
         open(3,file=fil1,status='unknown')

         print *,' Data file (MATLAB!!!!):'
         print *,fil1
         print *,'-------------------------------'
      endif


      call plotf(3,iform,t)
      close(3)
c ///////////////////////////////////////
c
      endif
         
c       if(mod(istep,nwrt).eq.0) then
         
!        if ((t-Tx).ge.tanim) then
!          Tx=t
!        iflag=0
c        if (mod(istep,nwrt*2).eq.0) iflag=1
!        call testdump(t,ybar,iflag,tstart)
       
!      endif

c        call stressb
c       endif

      twave = pi2/xomega

      if (t > waveperiods*twave) then
        write(*,*) 'Max. run duration of 10 wave periods exceeded !'
        stop
      endif

 1000 continue

c     timediff = dtime(timearray)
      write(*,*) timediff


      close(63)
      close(64)
      close(65)
      close(66)

      call timecall(isecf)
      write(*,*) 'TIMING FOR ENTIRE RUN:'
      call timeelapse(isecs,isecf)
      write(*,*) 'In secs.',
     >isecf-isecs
      call testhelm 
c     call testderv  

   15 format(a2,e10.4)
   25 format(a8,e10.4)
         stop
         end
