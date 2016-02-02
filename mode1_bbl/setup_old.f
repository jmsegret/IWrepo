C-This file contains all the routines related to setup of the code
C-I.e. wave-number definitions, filtering coefficients and Legendre
C-interpolation information.

c*********************************************************************  
      subroutine setup
c*********************************************************************  
      character*12 fin
      include 'parambctype.h'
      include 'dim.h'

c
      include 'comwave.h'
      include 'comleng.h'
      include 'comfilter.h'
      include 'comcoeff.h'
      include 'comzeta.h'
      include 'comwfunc.h'
      include 'comddzwfunc.h'
      include 'comlegendre.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comspecfilter.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'
      include 'comparamphys.h'
      include 'comtime.h'
      include 'comflow.h'
c     include 'comindex.h'
c                                                                       
c  determine the collocation and quadrature points in the z-direction
c  transform coordinates are used, i.e. z = [-1,1]
c

      parameter (niter=500)      

C-Set up collocation points      
      call jacobl(nzm,0.,0.,zpts,nzm)
C-Set up weights for use in weak formulation
C-Also set up legendre polynomials
      call quad(nzm,zpts,wg,nzm)
C-Set up derivative operator
      call derv(nzm,zpts,d,d2,d3,nzm)
C-Set up lengths and positions of subdomains
      call setup_subds

      open(206,file='rhomean.dat')
      write(206,*) 'VARIABLES = "`r_M(z)","z","N(z)","f(z)"'
C-Dump out stability curves for Robin BC penalty coefficient
      open(730,file='robin_penparam.dat')
      write(730,*) 'VARIABLES="`e","`t_m_i_n","`t_m_a_x"'
      do i=6,10
        fd = 10.**(float(i))
        epsilon = 1./fd
        ifbc = irobin
        ybar = zh(2)
c       call penalty_setup(epsilon)
c       call param_output(ifbc,alpha0,beta0,tau0,tau0max,
c    >                        epsilon,ybar,ks)
        write(730,*) epsilon,tau0,tau0max
      enddo
 
      close(730)

C-Calculate Collocation points on nz=65 grid
C-Use only when interpolating data from a finer (nz=97, e.g.)
C-grid
c     call jacobl(64,0.,0.,zpts2,64)
C-Calculate values of Legendre polynomials (based on Nz=97 interpolation)
C-on 64-point grid.
c     call quad2(64,zpts2,nzm)

c

c
c  setting the wave number
c     
c-x-dir. wavenumbers
      do 1 i=1,nxhp
         xw(i) = (i - 1)*alpha
         xsq(i) = xw(i)*xw(i) 
    1 continue

c-y-dir. wavenumbers
      do 2 j=1,nyh
         yw(j)       =   (j - 1)*beta 
         yw(ny+1-j)  =  -j*beta 
    2 continue
      do j=1,ny
         ysq(j) = yw(j)*yw(j)
      enddo

c     write(*,*) ysq(nyh+1),sqrt(xsq(nxhp)+ysq(nyh+1))
c     write(*,*) ysq(nyh),sqrt(xsq(nxhp)+ysq(nyh))
c     pause
c-z-dir. wavenumbers
c-(Only used for random noise generation)
      do 3 k=1,nzh
         zw(k)       =   (k - 1)*gamma
         zw(nzfl+1-k)  =  -k*gamma
    3 continue
      do k=1,nzfl
         zsq(k) = zw(k)*zw(k)
      enddo

c     write(*,*) sqrt((float(nxhp-1)*alpha)**2.+
c    +          (float(nyh-1)*beta)**2. +
c    +          (float(nzfl/2-1)*gamma)**2.),
c    + sqrt(xsq(nxhp)+ysq(nyh)+zsq(nzfl/2))

C-Calculate and Tabulate Boyd-Vandeven filter
      open(800,file='bvfilter.dat')
      write(800,*) 'VARIABLES ="MODE","BVFILTER","EXPFILTER"'
      kc =   fackc*nzm 
      kcbv = fackcbv*nzm
      pbv = 2.
      do k=0,nzm
        bvf = specfilter_bv(k,kcbv,nzm,pbv)
        write(800,'(1x,f10.5,2x,f15.10,2x,f15.10)') 
     >  float(k)/float(nzm),bvf,specfilter(k,kc,nzm,p,amp)
c       write(*,*) float(k)/float(nzm),specfilter(k,kc,nzm,p,amp)
      enddo
      close(800)

C-Calculate coordinate values of grid-points and
C-store in corresponding common block
      pi2 = 8.0*atan(1.0)
      xl = pi2/alpha
      yl = pi2/beta
      dy = yl/ny
      dx = xl/nx
c
      do i=1,nx
         x(i) = (i-1)*xl/nx
      enddo
c
      do j=1,ny
         y(j) = (j-1)*yl/ny
      enddo
c
C-Scan over subdomains (index ks)
C-Scan over all points in each subdomain (index kloc)
C-
      open(209,file='z_grid.dat',status='unknown')
      write(*,*) 'Origins and heights of each subdomain'
      do ks = 1,nsubd
c       write(*,*) 'Element start & thickness',z0(ks),zh(ks)
        do kloc = 1,nzloc
         zloc = 0.5*zh(ks)*(zpts(kloc-1)+1.0)
         k = (ks-1)*nzloc + kloc
         z(k) = z0(ks)+zloc
         write(209,*) k,z(k)
           
         k2 = (ks-1)*(nzloc-1) + kloc 
         zf(k2) = z0(ks)+zloc
!        write(*,*) k,z(k),k2,zf(k2),zpts(kloc-1)
c        write(*,*) k,zpts(kloc-1), 0.5*zh(ks)*(zpts(kloc-1)+1.0)
        enddo
        write(*,*) ks,z0(ks),zh(ks)  
      enddo
      close(209)


C-Calculate Mean Density Gradient
!      write(*,*) 'Mean Density Gradient Info Outputted'
! !     write(*,*) rhograd(z(nz),zlen,grav,rho0,brunt),brunt, 
!     >rhomean(zz,zlen,grav,rho0,brunt)
      do k=nz,1,-1
        zz = z(k)
!	  drhomean = abs( rhomean(z(nz),zlen,grav,rho0,brunt) -
!     >               rhomean(z(1),zlen,grav,rho0,brunt) )
!	  
        write(206,'(4(E12.6,2x))')
     >  rhomean(zz,zlen,grav,rho0,brunt),zz,
     >  bruntz(zz,zlen,brunt),0
c     >  fshear(zz,zdim)
      enddo
      close(206)

!-PD-2/15/08: No longer needed for critical layer problem     
!     write(*,*) 'Setting Up Wave Forcing' 
!     call wavesetup

c     goto 300        

c     
c-Set up coefficients used in filtering small wavenumbers
c-Subroutine taken from Senthil
c-PD-7/9/02: This is a slight modification of the global version.
c-We work with the global coordinates z(...) and those of the collocation
c-points ztps(...) in the [-1,1]. The spatial filtering doesn't care
c-about whether the global domain is partitioned into subelements or not.
C-NOTE: zpts indexing is: 0,...,nzm . z indexing is: 1,...,nz.
c
C-Careful: With DISCONTINUOUS PENALTY METHOD use zf(...)
C-coordinate

      do k=2,(nz2+1)/2-1
        kp=k+1
        km=k-1
        atpz(k)=(zf(k)-zf(km))/(zf(kp)-zf(km))/2.
        btpz(k)=( 2.+(2.*zf(kp)-4.*zf(k)+2.*zf(km))/
     >          (zf(kp)-zf(k)) )/4.
        ctpz(k)=2.*(zf(k)-zf(km))**2./
     >        (zf(kp)-zf(km))/(zf(kp)-zf(k))/4.
      end do

      do k=(nz2+1)/2+1,nz2-1
        kp=k+1
        km=k-1
        atpz(k)=2.*(zf(k)-zf(kp))**2./
     >        (zf(kp)-zf(km))/(zf(k)-zf(km))/4.
        btpz(k)=( 2.+(4.*zf(k)-2.*zf(kp)-2.*zf(km))/
     >          (zf(k)-zf(km)) )/4.
        ctpz(k)=(zf(kp)-zf(k))/(zf(kp)-zf(km))/2.
      end do

      atpz(1)=0.
      btpz(1)=1.
      ctpz(1)=0.

      atpz(nz2)=0.
      btpz(nz2)=1.
      ctpz(nz2)=0.

      atpz((nz2+1)/2)=1./4.
      btpz((nz2+1)/2)=2./4.
      ctpz((nz2+1)/2)=1./4.

c     do k=1,nz2
c       write(*,*) k,atpz(k),btpz(k),ctpz(k)
c     enddo

c     pause

 300  continue

c-The following is necessary in the solution of the Poisson eqn.
c-It is used via the patching scheme found in chptr. 22 of Boyd's
c-book.

c     write(*,*) '***************************************************'
c     write(*,*) 'Calculating Homogeneous ODE '
c     write(*,*) 'Solutions in each subdomain'
c     write(*,*) '***************************************************'
c     itemp = 0
c     write(*,*) 'For Surrogate Pressure PHI'
c     write(*,*) 'For zeta+-0,u & v'
c     ff = 78.956835208715
c     call calchomog(slhzetauvl,slhzetauvr,ff,ispecfilter,niter,
c    >itemp)
c     call calchomog(slhwl,slhwr,0.,ispecfilter,niter,itemp)
C *************** A.M.A. May,28th, 2008 ***************
              Zlf=zdim*sqrt(-log(0.1)/xbcrit)
c              do k=1,nz
c               if(z(k).ge.(zcen-Zlf)) then
                KLi=indexfi(zcen-Zlf)
cc                exit
c               endif
c               enddo
c                do k=1,nz
c               if(z(k).ge.(zcen+Zlf)) then
                KUi=indexfi(zcen+Zlf)
c                exit
c               endif
c               enddo
c               write(*,*)'KLi',KLi,KUi

      return                                                            
      end                                                               




c**********************************************************************
      subroutine setup_subds    
c**********************************************************************
C-This is where the size of starting point and length of each spectral
C-subdomain in the vertical direction are determined.
!-This routine assumes uniform height subdomains
      include 'dim.h'
      include 'comsubd.h'

      write(*,*) '-------------------------------'
      write(*,*) 'SUBDOMAIN LOCATIONS AND HEIGHTS'
      do k=1,nsubd
        zfrac = 1./float(nsubd)
        z0(k) = (k-1)*zlen*zfrac
        zh(k) = zlen*zfrac
        isubdtype(k) = 1
        write(*,*) z0(k),zh(k)
      enddo

      write(*,*) '-------------------------------'
 

C-Information regarding types of subdomains of different height
      subdlen(1) = zlen*zfrac


      return
      end


!**********************************************************************
      subroutine setup_subd
!**********************************************************************
!-This is where the size of starting point and length of each spectral
!-subdomain in the vertical direction are determined.
      include 'dim.h'
      include 'comsubd.h' 
      include 'comcritlayer.h' 
      namelist/inputs/ istart,dtin,nsteps,nwrt,iform,cdtin
      namelist/params/ xlen,ylen,zlen,brunt,xnu,xkappa,grav,
     >expan,umax,rho0,fcshear,xacrit,xbcrit,zcen,umax

c      namelist/mflags/ irelax,ispecfilternl,
c     >ispecfilterexp,nspecfilterexp,amp,p,pf,fackc,fackcbv,
c     >facdir,facneumn,facrobin,navgintrfc,ncalclegspec,
c     >ifilter,nstepfilter,
c     >tstart,npostp,itempreset,nrestartdump,n3doutput,i0
c        real   Rz,Rz2,zrat
c        integer Nlo,Nup
      parameter(zmin= 0.06)
      logical c
      real rs(1000),remain(1000),remainder,rem0,rem1,remi,
     >     b0,b1,zmin,zlen,zsh,slope,b,icen,one,
     >     two,ratio



!
!-NOTE: These subdomains are referenced in a global (over the
!-      the entire computational domain) sense. Most other relevant
!-      references are done in a local sense (over a parallelization
!-      subdomain) (PD: 5/19/05).
!

!-PD-9/28/05: Subdomain information is now read in from a file

!-------------------------------------
!-Vertical *Computational* Subdomains
!- Calculate Computational subdomains instead of reading in the file "subd.dat"
!A.M.A. Feb-21-2008 when h is uniform (uniform subdomain heights)
!Note nsubd > 1
!-------------------------------------
c        open(987,file='testsetyup')
c             write(987,*)zcen

       if(nsdtype.eq.1) then    
           nsubd_in=nsubd
           z0(1)=0.0
          zh(1)=zlen/nsubd
          do ks=1,nsubd-1
            zh(ks+1)=zh(1)
            z0(ks+1)=z0(ks)+zh(1)
          enddo 
c***********This part of code clusters the grid around the localization
c function such that the grid size is minimum at Zsh and expands in both
c directions, away from Zsh, with expansion ratio Rz
c A.M.A. May,21,2008
             elseif(nsdtype.eq.2)then
c           open(700,file='subd.dat',status='old')
c           read(700,*) nsubd_in
c           write(*,*) 'Number of subdomains',nsubd_in
c           do ks=1,nsubd
c            read(700,*) z0(ks),zh(ks)
c           enddo    
c           close(700)
c c           endif

! remove first column of comments for Ammar's subdomain spacer

CC&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c               zsh2=zsh+0.1
c               Nlo=5 
c               Nup=nsubd-Nlo              
c               Rz=4.0**(1.0/(Nup-1))      
c              zh(Nlo+1)=(zlen-zsh2)*(Rz-1.0)/((Rz**Nup)-1.0)
c              z0(Nlo+1)=zsh2
ccC*************Setting size of sudomains above the c.l.*********
c               do  k=2,Nup
c               zh(k+Nlo)= zh(k+Nlo-1)*Rz
c               z0(k+Nlo)=z0(k+Nlo-1)+zh(k+Nlo-1)
c               enddo             
cC*************Setting size of subdomains below the c.l.********
cc               zh(Nlo)=1.3*zh(Nlo+1)
cc               z0(Nlo)=z0(Nlo+1)-zh(Nlo)
cC-------------------------------------------- 
cc               zh(Nlo-1)=1.3*zh(Nlo)
cc               z0(Nlo-1)=z0(Nlo)-zh(Nlo-1)
cC-----------------------------------------------
cc               zh(Nlo-2)=1.3*zh(Nlo-1)
cc               z0(Nlo-2)=z0(Nlo-1)-zh(Nlo-2)
cC----------------------------------------------------
cc               zh(2)=0.5*(zsh2-(zh(Nlo)+zh(Nlo-1)+zh(Nlo-2)))
cc               z0(2)=z0(3)-zh(2)
cc               zh(1)=zh(2)
cc               z0(1)=0.0
cc
c               zrat=4  !(zh(nsubd)*zsh2)/(zh(Nlo+1)*(zlen-zsh2))
c               Rz2=zrat**(1.0/(Nlo-1))
c               zh(Nlo)=zsh2*(Rz2-1.0)/((Rz2**Nlo)-1.0)
c               z0(Nlo)=zsh2-zh(Nlo)
c                z0(1)=0.0
c                zh(1)=zh(Nlo)*zrat
c              do k=2,Nlo-1
c                z0(k)=z0(k-1)+zh(k-1)
c                zh(k)=zh(k-1)/Rz2
c              enddo
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

C MPR Jun 2010
C This part localizes the vertical subdomain around zsh
C with each subdomain being some constant ratio of the one
C before it. The start of each subdomain follows
C a geometric series decreasing from z = 0 to zsh, and follows
C the same geometric series increasing from zsh to zlen. 
C zmin sets the first subdomain thickness. ADJUST IT ABOVE! ^^

      write(*,*) 'SOLVING SUBDOMAIN'


      if((zmin.lt.(zlen/nsubd)).or.(zmin.gt.(0.4*zsh))) then    
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*****SUBDOMAIN INPUTS NOT VALID******'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      write(*,*) '*************************************'
      endif

C This part takes a snapshot of the function we want to zero.
C The function, 'remainder', calculates how close the final
C subdomain location is to zlen, given some ratio of
C subdomains, held in the vector 'rs'


      open(10,file='rem.dat')
      open(12,file='rs.dat')
      do i = 1,1000
          rs(i) = 0.3 + i*0.00069
          write(12,*) rs(i) 
      enddo
      close(12)

      do i = 1,1000
          remain(i) =
     >  remainder(rs(i),(10.0*zmin/zlen),nsubd,10.0,(10.0*zsh/zlen))
          write(10,*) rs(i), remain(i)
c, abs(remain(i))
      enddo

      close(10)

C Since 'remainder' has multiple roots, we want to find the 
C highest valued root. The largest root is the only root that is 
C not complex. This function searches 'remainder' backwards
C until it hits a zero. This first estimate is where
C we start the Newton-Raphson method to make sure
C we get the highest valued root.

      rem1 = 1.0e20
      c = 0
      i = 1000

c      open(20,file='fguess.dat')

      do while(c.eq.0)
           remi = remain(i)
           rem0 = abs(remi)
           if((rem0.lt.rem1).or.(rem0.gt.7)) then
                spot = i
                c = 0
                rem1 = rem0
c                write(20,*) rem1, rem0
           elseif(rem0.gt.rem1) then
                c = 1
           endif
           i = i - 1
      enddo

      if((c.eq.0).or.(i.gt.997)) then
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****SUBDOMAIN LOCALIZATION FAILURE******'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      write(*,*) '*****************************************'
      endif

c      write(20,*) '*************************'
c      write(20,*) spot, rs(spot)
c      write(20,*) '*************************'
c      close(20)

C This is part uses Newton-Raphson to get a good
C close estimate of the ratio of subdomains, b.

      b0 = rs(spot) + 0.0001
      b1 = b0 + 0.00001

      rem0 = remainder(b0,zmin,nsubd,zlen,zsh)
      rem1 = remainder(b1,zmin,nsubd,zlen,zsh)

      do while(abs(rem1).gt.0.000001)

           slope = (rem1 - rem0)/(b1 -b0)
           b0 = b1
           b1 = b1 - (rem1/slope)

           rem0 = rem1

           rem1 = remainder(b1,zmin,nsubd,zlen,zsh)
      enddo

      b = b1

      icen = log(abs(1. - (zsh*(1.-b)/zmin)))/(log(b))

C This part calculates subdomain locations

c      open(630,file='z0.dat')
      z0(1) = 0.0
      do i = 1,(nsubd-1)
           if(i.lt.icen) then
                z0(i+1) = zmin*(1. - (b**i))/(1. - b)
c                write(630,*) i, z0(i)
           else
                one = 2.*(1. - b**(icen))/(1. - b)
                two = (1. - b**(2.*icen - i))/(1. - b)
                z0(i+1) = zmin*(one - two)
c                write(630,*) i, z0(i)
           endif
      enddo
c      write(630,*) nsubd,z0(nsubd)

c      close(630)

      ratio = z0(2)/(z0(2)-z0(3))

      write(*,*) 'Subdomain ratio is',ratio

c*********************** write to a file ***********
         open(700,file='subd.dat')
         write(700,*)nsubd
            do ks=1,(nsubd-1)
             zh(ks) = (z0(ks+1) - z0(ks))
             write(700,*) z0(ks),zh(ks)
            enddo
            zh(nsubd) = (zlen - z0(nsubd))
            write(700,*) z0(nsubd),zh(nsubd)
         close(700)

c          write(*,*)'zsh=',zsh

c************END Of Grid Generation ********************
          else
          open(700,file='subd.dat',status='old')
          read(700,*) nsubd_in
          write(*,*) 'Number of subdomains',nsubd_in
          do ks=1,nsubd
           read(700,*) z0(ks),zh(ks)
          enddo    
          close(700)
         endif

      end subroutine  setup_subd










c*******************************************************
c Remainder is used in setup_subd if nsdtype is not 1

      function remainder(r,zmin,nsubd,zlen,zsh)
      implicit none
      integer nsubd
      real r,remainder,zmin,zlen,zsh,a,icen,one,two


      a = zmin
      icen = log(abs(1. - (zsh*(1. - r)/a)))/(log(r))
      one = 2.*(1. - r**(icen))/(1. - r)
      two = (1. - r**(2.*icen - nsubd))/(1. - r)
      remainder = a*(one - two) - zlen

      end



c**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
c**********************************************************************
C----------------------------------------------------------------------
C-Calculates values of all Legendre polynomials (and immediately highest
C-one in hierarchy) at a given collocation point xc.
C-Same operation for derivatives of Legendre polynomials.
      dimension al(0:ndim),alp(0:ndim)
c
      al(0) = 1.
      al(1) = xc
      alp(0) = 0
      alp(1) = 1.
c
      do k=1,n
        kp = k + 1
        km = k - 1
        al(kp) = (2.*k+1.)*xc*al(k)/kp - k*al(km)/kp
      enddo
c
      do k=1,n
        kp = k + 1
        km = k - 1
        alp(kp) = (2.*k+1.)*(al(k)+xc*alp(k))/kp - 
     &            k*alp(km)/kp
      enddo
c
      return
      end
c

c**********************************************************************
      subroutine legenaux(alaux,n,xc,ndim)
c**********************************************************************
C----------------------------------------------------------------------
C-Special addition by PD on 4/5/04: Routine that calculates values
C-of single domain Legendre polynomials on multidomain grid.
C-(Used to interpolate TG-eqn. eigenfunction on multidomain grid).
C-
C-Calculates values of all Legendre polynomials (and immediately highest
C-one in hierarchy) at a given collocation point xc.
C-Same operation for derivatives of Legendre polynomials.
      dimension alaux(0:ndim)
c
      alaux(0) = 1.
      alaux(1) = xc
c
      do k=1,n
        kp = k + 1
        km = k - 1
        alaux(kp) = (2.*k+1.)*xc*alaux(k)/kp - k*alaux(km)/kp
      enddo
c
      return
      end
c





c**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
c**********************************************************************
c
c  computes the gauss-lobatto collocation points
c
c   N:              degree of approximation (order of polynomials=n+1)
c   ALPHA:          parameter in Jacobi weight
c   BETA:           parameter in Jacobi weight
c   XJAC:           roots from largest to smallest
c
c
c  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
c  for Legendre-Gauss-Lobatto points use           0             0
c
      dimension xjac(1000),xcol(0:ndim)
      common /jacpar/alp,bet,rv
      data kstop/10/
      data eps/1.0e-12/
c
      alp = alpha
      bet = beta
      rv = 1. + alp
      np = n + 1
c
        
      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1.0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1.0)
      det = pnp*pnm1m - pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m
      a = (rp*pnm1m - rm*pnm1p)/det
      b = (rm*pnp   - rp*pnm)/det
      xjac(1) = 1.0
      nh = (n+1)/2
      dth = 3.14159265/(2*n+1)
      cd = cos(2.*dth)
      sd = sin(2.*dth)
      cs = cos(dth)
      ss = sin(dth)
c
      do 39 j=2,nh
       x = cs
       do 29 k=1,kstop
        call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1 + a*pn + b*pnm1
        pder = pdnp1 + a*pdn + b*pdnm1
        recsum = 0.0
        jm = j - 1
        do 27 i=1,jm
          recsum = recsum + 1.0/(x-xjac(i))
27      continue
28      continue
        delx = -poly/(pder-recsum*poly)
        x = x +delx
        if(abs(delx) .lt. eps) go to 30
29      continue
30      continue
        xjac(j) = x
        cssave = cs*cd - ss*sd
        ss = cs*sd + ss*cd
        cs = cssave
39      continue
        xjac(np) = -1.0
        npp = n + 2
        do 49 i=2,nh
          xjac(npp-i) = -xjac(i)
49      continue
        if(n.ne. 2*(n/2)) then 
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         go to 56
        else
        xjac(nh+1) = 0.0
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         return
        endif

c
56      return
        end
c




c**********************************************************************
        subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
c**********************************************************************
        common /jacpar/ alp,bet,rv
        apb = alp + bet
        poly = 1.0
        pder = 0.0
        if(n.eq.0) return
        polylst = poly
        pderlst = pder
        poly = rv*x
        pder = rv
        if(n.eq.1) return
        do 19 k=2,n
          a1 = 2.*k*(k+apb)*(2.*k+apb-2.)
          a2 = (2.*k+apb-1.)*(alp**2-bet**2)
          b3 = (2.*k+apb-2.)
          a3 = b3*(b3+1.)*(b3+2.)
          a4 = 2.*(K+alp-1.)*(k+bet-1.)*(2.*k+apb)
          polyn = ((a2+a3*x)*poly - a4*polylst) / a1
          pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
          psave = polylst
          pdsave = pderlst
          polylst = poly
          poly = polyn
          pderlst = pder
          pder = pdern
19      continue
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
        return
        end
c






c**********************************************************************
      subroutine quad(n,x,w,ndim)
c**********************************************************************
      include 'dim.h'
      include 'comlegendre.h'
      include 'comsubd.h'
      parameter (nn=1000)
      dimension x(0:ndim),w(0:ndim),alp1(0:nn),al1(0:nn)
      dimension xaux(0:nzaux-1),alaux1(0:nn)
c
c  **PD-modify**: Also store values of 0 to ndim Legendre
c  polynomials at all collocation points.
c  determine the Gauss Quadrature weighting factors
c
      small = 1.0e-30
      do k=0,n
       xc = x(k)
       call legen(al1,alp1,n,xc,nn)  
c-Calculate weighting factor
       w(k) = 2. / 
     &         ( n*(n+1)*al1(n)*al1(n) + small )
c-Store values for j-th Legendre polynomial at point xc=x(k)
c-Store values for derivative of j-th order L. polynomial at
c-collocation point, as well.
       do j=0,nzm
         alpol(j,k) = al1(j)
         dalpol(j,k) = alp1(j)
       enddo

      enddo

c     open(210,file='legen_out.dat')
c     do k=1,nz
c       write(210,*) x(k-1),alpol(16,k-1)
c     enddo
c     close(210)

      goto 100

C-Special addition: Reads in multidomain grid and interpolates
C-values of Legendre polynomials on this grid
C-Read in file with multidomain coordinates
C-USED ONLY IN SINGLE DOMAIN CALCULATIONS
      open(209,file='z_grid.dat',status='unknown')
      do k=1,nzaux
        read(209,*) idummy,xin
        zaux(k) = xin*.25
        xaux(k-1) = (2.*xin-zlen)/zlen
c       write(*,*) idummy,xin,xaux(k-1)
      enddo

C-Now interpolate values of single domain Legendre polynomials
C-PD: After introducing simpler interpolation scheme
C-leave this for future use (4/6/04)
      do k=0,nzaux-1
c       write(*,*) k,xaux(k)
        xc = xaux(k)
        call legenaux(alaux1,n,xc,ndim)

        do j=0,nzm
          alpolaux(j,k) = alaux1(j)
        enddo
 
      enddo

c-savva
c     open(210,file='legen_out.dat')
c     do k=1,nzaux
c       write(210,*) xaux(k-1),alpolaux(16,k-1)
c     enddo
c     close(210)
      
  100 continue
      return
      end

c**********************************************************************
      subroutine quad2(n,x,ndim)
c**********************************************************************
C-USED ONLY WHEN WANTING TO PROJECT TO LOWER VERTICAL RESOLUTION FIELD
      include 'dim.h'
      include 'comlegendre.h'
      parameter (nn=1000)
      dimension x(0:ndim),alp1(0:nn),al1(0:nn)
c
c  **PD-modify**: Also store values of 0 to ndim Legendre
c  polynomials at all collocation points (for Nz=65 grid).
c
C-This routine is only useful when trying to read in a field interpolated
C-from a more fine grid-direction
      small = 1.0e-30
      do k=0,n
       xc = x(k)
       call legen(al1,alp1,n,xc,nn)
c-Store values for j-th Legendre polynomial at point xc=x(k)
       do j=0,nzm
         alpol2(j,k) = al1(j)
       enddo

      enddo

      return
      end
c





c**********************************************************************
      subroutine derv(nterm,x,d,d2,d3,ndim)
c**********************************************************************
      parameter (nn=1000)
      dimension x(0:ndim),
     >          d(0:ndim,0:ndim),d2(0:ndim,0:ndim),d3(0:ndim,0:ndim),
     >          al1(0:nn),alp1(0:nn),
     >          al2(0:nn),alp2(0:nn),c(0:nn)
c
c  determine the derivative at the collocation points
c  The following is the conventional approach
c  which generates errors, particularly at higher
c  order derivative evaluation (see Costa and Don).

C-**PD: 2/14/03**:
C-USE THIS ONLY FOR CALCULATION OF LEGENDRE POLYNOMIALS
C-AT COLLOCATION POINTS
      do i=0,nterm
        xi = x(i)
        call legen(al1,alp1,nterm,xi,nn)  
       do j=0,nterm
        xj = x(j)
        call legen(al2,alp2,nterm,xj,nn)  
c       if(i.eq.j) then
c        d(i,j) = 0
c       else
c        d(i,j) = al1(nterm)/(al2(nterm)*(xi-xj))
c       endif
       enddo
      enddo

c
c     ann = 0.25*nterm*(nterm+1)
c     d(0,0) = -ann
c     d(nterm,nterm) =  ann

C-Introduce methodology by Costa and Don.
C-Calculate off diagonal elements of D^1
C-1) Must first compute coefficients C_i
      do 10 k=0,nterm
        prod = 1.
        xk = x(k)

        do 20 l=0,nterm
          xl = x(l)

          if (l.ne.k) then
            prod = prod*(xk-xl)
          else
            prod = prod
          endif

 20     continue
        c(k) = prod

 10   continue

C-This is the lengthy route, but no big deal.
C-FIRST DERIVATIVE   
C-Calculate off diagonal elements.
      do 30 k=0,nterm
        xk = x(k) 

        do 40 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d(k,j) = c(k)/(c(j)*(xk-xj))
          else
            d(k,j) = 0.
          endif

 40     continue
 30   continue

C-Calculate diagonal elements now
C-Diagonal element is negative row sum of off diagonal
C-elements (See Costa & Don)      
      do 50 k=0,nterm
        sum = 0.
        do 60 j=0,nterm

          if (k.ne.j) then
            sum = sum + d(k,j)
          else
            sum = sum
          endif

 60     continue
        d(k,k) = -sum
c       if ((k.ne.0).and.(k.ne.nterm)) d(k,k) = 0.
 50   continue

C-SECOND  DERIVATIVE
C-Off-diagonal elements
      m = 2
      do 70 k=0,nterm
        xk = x(k)

        do 80 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d2(k,j) = float(m)*(d(k,k)*d(k,j) - d(k,j)/(xk-xj))
          else
            d2(k,j) = 0.
          endif

 80     continue
 70   continue      

C-Diagonal elements
       do 90 k=0,nterm
        sum = 0.
        do 100 j=0,nterm

          if (k.ne.j) then
            sum = sum + d2(k,j)
          else
            sum = sum
          endif

 100    continue
        d2(k,k) = -sum
 90   continue

C-THIRD  DERIVATIVE
C-Off-diagonal elements
      m = 3
      do 110 k=0,nterm
        xk = x(k)

        do 120 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d3(k,j) = float(m)*(d2(k,k)*d(k,j) - d2(k,j)/(xk-xj))
          else
            d3(k,j) = 0.
          endif

 120    continue
 110  continue

C-Diagonal elements
       do 130 k=0,nterm
        sum = 0.
        do 140 j=0,nterm

          if (k.ne.j) then
            sum = sum + d3(k,j)
          else
            sum = sum
          endif

 140    continue
        d3(k,k) = -sum
 130  continue

c     do i=0,nterm
c       do j=0,nterm
c         write(*,*) i,j,d(i,j),d2(i,j),d3(i,j)
c       enddo
c     enddo

      return
      end
c

c******************************************************************
      function zderivmap(ybar)
c******************************************************************
C
C-PD: 2/21/02. Real function that returns value of z-derivative
C-(dzeta(z)/dz) of vertical mapping function which maps a domain
C-of the type (z0,z0+ybar) to the domain (-1,1) used in Legendre interpolation.
C-For now the derivative of the mapping function is simply a constant
C-(See notes)
C
      zderivmap=2./ybar
 
      return
      end
 
 
c******************************************************************
      function zderivmap2(ybar)
c******************************************************************
C
C-Coefficient for calculating 2nd derivative
C-(See notes)
C
      zderivmap2=(2./ybar)**2.

      return
      end


c******************************************************************
      function zderivmap3(ybar)
c******************************************************************
C
C-Coefficient for calculating 2nd derivative
C-(See notes)
C
      zderivmap3=(2./ybar)**3.

      return
      end
