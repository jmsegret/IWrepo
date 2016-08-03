c****************************************************************************  
      subroutine init ! (zlen)
c*********************************************************************
      include 'dim.h'
      parameter (initspec=0,np=25)
c     
      include 'comflow.h'
      include 'comwave.h'
      include 'comtime.h'
      include 'comleng.h'
      include 'combc.h'
      include 'comscratch.h'
      include 'comparamphys.h'
      include 'comtkeprev.h'
      include 'comgridptxyz.h'
      include 'comrms.h'
      include 'comsubd.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'

      dimension unoise(nxpp,ny,nz)
      dimension vnoise(nxpp,ny,nz)
      dimension wnoise(nxpp,ny,nz)
c     +temp(nxpp,ny,nz)
      dimension ozz(nxpp,ny,nz)
      dimension spec(3,np),spectmp(np)
c     dimension utemp(nz),vtemp(nz),wtemp(nz)
c     dimension x(nx),y(ny),z(nz) ! ,u1(ntot),z2(nz)
c     dimension finu(nzloc),foutu(nzloc) 
c     dimension finv(nzloc),foutv(nzloc)
c     dimension finw(nzloc),foutw(nzloc)
      complex rx 

c///////////////////////////////////////////
      
c     dimension  x(nx),y(ny),z(nz)
      rx=-1.0*cmplx(0.0,1.0)

      open(552,file='noisespectrum.dat')
      open(551,file='rmsprofilesz.dat')
      open(550,file='rmsprofilesy.dat')

c
c setup the vertical uniform mesh 
c      
c     

c     call testftv

C-Set Initial Condition (Absolutely no perturbations
C-=PD:2/25/04)

c       xomegandim = (xomega*zdim)/(umax)

      do i=1,nxpp
        do j=1,ny
          do k=1,nz
c       sinf = sin( xkcrit*x(i) + xmcrit*z(k))
c       cosf = cos( xkcrit*x(i) + xmcrit*z(k))
          
c       undim=(-(xacrit*xmcrit/xkcrit)*fzloc(z(k),zcen,zdim,xbcrit)
c     >    *cosf)-((xacrit/xkcrit)*fzlocp(z(k),zcen,zdim,xbcrit)*sinf)       

c       u(i,j,k)=umax*undim
 
           u(i,j,k) = 0.                                                
           v(i,j,k) = 0.                                                  

c           wndim = xacrit*fzloc(z(k),zcen,zdim,xbcrit)*cosf
c           w(i,j,k)=umax*wndim

           w(i,j,k) = 0.                                                  
c        xomegandim = (xomega*zdim)/(umax)
     
c        tempndim = -(xacrit/xomegandim)*fzloc(z(k),zcen,zdim,xbcrit)
c     >               *sinf
c        temp(i,j,k)=tempndim*zdim*brunt*brunt*rho0/grav
      
           temp(i,j,k) = 0.                                               


           ox(i,j,k) = 0.
           oy(i,j,k) = 0.
           oz(i,j,k) = 0.
           ozz(i,j,k) = 0.
          
          enddo
        enddo
      enddo

      do i=1,np
        spectmp(i) = 0. 
      enddo

c     
c---------------------------------------------------------------
c     -
c     linear temperature profile                                  -
c     in the physical space                                       -
c     -
c---------------------------------------------------------------
c     
c     imposed with random fluctuations
c     
c     for Cray
      
c     call ranset(57005)
c     
      do k=1,nz
         do j=1,ny       
            do i=1,nx

c-PD: Linear Stratification 
C-PD: 5/28/03: There has been a slight change here. We
C-are now working with the temperature PERTURBATION
               zz = z(k)
               temp(i,j,k)=  0. 

c              if ((i.eq.nxh).and.(j.eq.nyh)) write(*,*) k,zz,temp(i,j,k)
            enddo
         enddo
      enddo
c     /////////////////////////////////////////////

      pi2 = 8.0*atan(1.0)
      xl = pi2/alpha
      yl = pi2/beta
!      write(*,*) 'umax',umax      

      do k=1,nz
        do j=1,ny
          do i=1,nx
     
            ox(i,j,k) = 0.
            oy(i,j,k) = 0.
            oz(i,j,k) = 0.
 
          enddo
        enddo
      enddo

c////////////////////////////////////////////////
c
c
c     skip setup for flux boundary conditions
c     go to 99
c
c     impose temperature gradient (= -heat flux/kappa)
c     along upper and lower walls
c     Flux setup here !!!!

      do j=1,ny
         do i=1,nx
      
c     /////// convection setup
c     
c     ox(ijk) = -1.0
c     ox(ij)  = -1.0 
            
c           ox(ijk) = 0.0
c           ox(ij)  = 0.0 

c-Constant flux temperature BC's for stratified flow
!           ox(i,j,nz) = 1.0
!           ox(i,j,1) = 1.0
c-Fixed temperatures at boundaries
            ox(i,j,nz) = 0.0
            temp(i,j,nz)=0.0 
            ox(i,j,1)  = 0.  
            temp(i,j,1)= 0.  

c-Set BC values for velocity field (oy-u, oz-v, ozz-w)
c-Assumes solid bottom wall and free slip upper surface.
c-u-perturbation at bottom is such that soliton should
c-match free stream (remember, moving reference frame)

! JMS-edit try to get bc correct

            oy(i,j,1) = 1.0
            oy(i,j,nz) = 1.0
            oz(i,j,1) = 0.
            oz(i,j,nz) = 0.
            ozz(i,j,1) = 0.
            ozz(i,j,nz) =  0. 

         enddo
      enddo
      
      
c--------------------------------------------------
c     transform to spectral space                                          
c--------------------------------------------------
c     
     

      call  horfft (u,-1,nz)
      call  horfft (v,-1,nz)
      call  horfft (w,-1,nz)

     
c----------------------------------------------------------------------
      
c     
      
      call  horfft (temp,-1,nz)
      call  horfft (ox,-1,nz)

      call  horfft (oy,-1,nz)
      call  horfft (oz,-1,nz)
      call  horfft (ozz,-1,nz)

c     
c     save the boundary condition along the upper boundary (tubc)
c     and along the bottom boundary (tbbc) (for temperature)

c     Apply same process for u,v & w variables
c

C-**PD**2/18/02: WHY NOT DO THIS DIRECTLY ? WHY USE ox ARRAY ?
C-It is necessary when one needs to convert BC's to spectral space
      do j=1,ny
         do i=1,nxpp
            tubc(i,j) = ox(i,j,nz)
            tbbc(i,j) = ox(i,j,1)

            uw(i,j,2) = oy(i,j,nz)
            uw(i,j,1) = oy(i,j,1)

            vw(i,j,2) = oz(i,j,nz)
            vw(i,j,1) = oz(i,j,1)

            ww(i,j,2) = ozz(i,j,nz)
            ww(i,j,1) = ozz(i,j,1)
         enddo
      enddo
c     
 99   continue             

      close(552)
      close(551)
      close(550)
     
      return                                                            
      end                                                               




C******************************************************************* 
      subroutine noise(fnoise3d,spec)
C*******************************************************************
C-Created by PD: 4/08/03.
C-The purpose of this routine is to generate a field in 3-D Fourier
C-space which is essentially a k^-5/3 spectrum with a slight perturbation
C-to it. The z-direction is considered periodic strictly for the
C-purpose of noise generation.

      include 'dim.h'

C-ATTENTION: Cut-off and max. wavenumbers are based assuming same number of
C-of Fourie rmodes in each direction
      parameter (kc=8,kmaxx=nxhp,kmaxy=nyh,kmaxz=nzh,kmaxl=nzh)
 
      parameter (expk=-11./6.,amp=10.) 
      parameter (np=25)

      include 'comwave.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comscratch.h'

      dimension fnoise3d(nxpp,ny,nz)

      dimension foutr(nz),fouti(nz)
      dimension spec(np),nspec(np)
      complex sum,ampf,expf,rx,tmp,cmode
      rx=-1.0*cmplx(0.0,1.0) 

      open(555,file='noiseout.dat')

      pi2 = 8.0*atan(1.0)

C-Set cut-off wavenumber and maximum wavenumber magnitudes
      akc = sqrt(zsq(kc))
c     akmax = sqrt(3.)*sqrt(xsq(kmax))
      akmax = sqrt(zsq(kmaxz))
      akmax2 = sqrt(xsq(nxhp) + ysq(nyh) + zsq(nzfl/2))
      akmaxp = 300.
c     write(*,*) 'Max wave #',akmax,akmax2
       

C-Nulify spectra
      do i=1,np

        spec(i) = 0.
        nspec(i) = 0
 
      enddo

c-Generate noise
      ii = 0
      do 100 i=1,nxpp,2
         ii = ii+1
         do 200 j=1,nyh
           jneg = ny+1-j

c          write(*,*) ii,j
C-Scan over all Fourier ghost points (these are simply needed to then
C-project on the Lobatto-Legendre grid)
           do 300 k=1,nzh
             kneg = nzfl+1-k

C-NOTE: We need to cover all 4 quadrants on the Okykz plane in
C-wavenumber space. (ak_i = indicates wavenumber magnitude in
C-corresponding quadrant)
         
             ak1 = sqrt(zsq(k))
             ak2 = sqrt(zsq(k))
             ak3 = sqrt(zsq(kneg))
             ak4 = sqrt(zsq(kneg))

C-Generate random perturbation coefficients for each velocity
C-component

c     if ((ak1.gt.akc).and.(ak1.le.akmax)) write(*,*) ii,j,k,ak1,akmax
c     if ((ak2.gt.akc).and.(ak2.le.akmax)) write(*,*) ii,j,k,ak2,akmax
c     if ((ak3.gt.akc).and.(ak3.le.akmax)) write(*,*) ii,j,k,ak3,akmax
c     if ((ak4.gt.akc).and.(ak4.le.akmax)) write(*,*) ii,j,k,ak4,akmax
     
C-For spectrum computation

C-Create noise functions
c            if ((ii.eq.1).and.(j.eq.1)) write(*,*) cr,ci

C-In 1st quadrant on Okykz plane

c            write(*,*) 'BOB',amp
             ampmode = ampnoise(ak1,akc,akmax,expk,amp)
             arg = random()*pi2*0.5
             phaser = cos(arg)
             phasei = sin(arg)
             cmode = ampmode*cmplx(phaser,phasei)
             fnoise3d(i,j,k) = real(cmode)
             fnoise3d(i+1,j,k) = real(rx*cmode)

c            nbin = int(ak1/akmaxp*float(np))+1
c            nspec(nbin) = nspec(nbin) + 1
c            spec(nbin) = spec(nbin) + cabs(cmode)**2.

C-In 2nd quadrant on Okykz plane

             ampmode = ampnoise(ak2,akc,akmax,expk,amp)
             arg = random()*pi2*0.5
             phaser = cos(arg)
             phasei = sin(arg)
             cmode = ampmode*cmplx(phaser,phasei)
             fnoise3d(i,jneg,k) = real(cmode)
             fnoise3d(i+1,jneg,k) = real(rx*cmode)

c            nbin = int(ak2/akmaxp*float(np))+1
c            nspec(nbin) = nspec(nbin) + 1
c            spec(nbin) = spec(nbin) + cabs(cmode)**2.

C-The following kz-wavenumbes are the highest ones and
C-cause trouble. Drop them.
C-In 3rd quadrant on Okykz plane
             ampmode = ampnoise(ak3,akc,akmax,expk,amp)
             arg = random()*pi2*0.5
             phaser = cos(arg)
             phasei = sin(arg)
             cmode = ampmode*cmplx(phaser,phasei)
             fnoise3d(i,jneg,kneg) = real(cmode)
             fnoise3d(i+1,jneg,kneg) = real(rx*cmode)

c            nbin = int(ak3/akmaxp*float(np))+1
c            nspec(nbin) = nspec(nbin) + 1
c            spec(nbin) = spec(nbin) + cabs(cmode)**2.

C-In 4th quadrant on Okykz plane
             ampmode = ampnoise(ak4,akc,akmax,expk,amp)
             arg = random()*pi2*0.5
             phaser = cos(arg)
             phasei = sin(arg)
             cmode = ampmode*cmplx(phaser,phasei)
             fnoise3d(i,j,kneg) = real(cmode)
             fnoise3d(i+1,j,kneg) = real(rx*cmode)

c            nbin = int(ak4/akmaxp*float(np))+1
c            nspec(nbin) = nspec(nbin) + 1
c            spec(nbin) = spec(nbin) + cabs(cmode)**2.


 300       continue
 200     continue
 100  continue

C-Now scan over all Lobatto Legendre points, defined in set-up.
C-Use eqn. (2.1.24) of Canuto et al. to project vertical
C-Fourier series on  LL points.
C-k=Index for LL point
C-kf = index for vertical Fourier mode
      ii = 0
      do 150 i=1,nxpp,2
         ii = ii+1
         do 250 j=1,ny

           do 350 k=1,nz

C-At each (i,j,k) point in horizontal Fourier space and
C-vertical physical space calculate noise content
             sum = cmplx(0.,0.)

             do kf=1,nzh
               kneg = nzfl+1-kf

               ak1 = sqrt(xsq(ii) + ysq(j) + zsq(kf))
c              ak2 = sqrt(xsq(ii) + ysq(jneg) + zsq(k))
c              ak3 = sqrt(xsq(ii) + ysq(jneg) + zsq(kneg))
c              ak4 = sqrt(xsq(ii) + ysq(j) + zsq(kneg))

               zz = (z(k)/zlen)*pi2
C-Contribution from positive kz-wavenumbers
               if (kf.lt.kmaxl) then
C-Careful about how you define the wavenumber
C-Argument is: index*2pi*(k/N) (i.e. the position of the grid-point
C-over an interval [0,1]).
                  arg = (kf-1)*zz               
                  expf = cmplx( cos(arg),sin(arg) )
                  ampf = 
     >            cmplx( fnoise3d(i,j,kf),fnoise3d(i+1,j,kf) )
                  sum = sum + ampf*expf 
C-Contribution from negative kz-wavenumbers
c              if ((kneg).lt.kmaxl) then
                  arg = -kf*zz
                  expf = cmplx( cos(arg),sin(arg) )
                  ampf = 
     >            cmplx( fnoise3d(i,j,kneg),fnoise3d(i+1,j,kneg) )
                  sum = sum + ampf*expf
               else
                  ampf = 0.
                  expf = cmplx(0.,0.)
               endif
               
             enddo 

             foutr(k) = real(sum)
             fouti(k) = real(rx*sum)

 350       continue
  
           do k=1,nz
            
             fnoise3d(i,j,k) = foutr(k)
             fnoise3d(i,j,k+1) = fouti(k) 
 
           enddo


 250     continue
 150  continue         

C-Transform noise fields to physical space
      call  horfft (fnoise3d,1,nz)

C-Dump noise Fourier spectrum file
c     do i=1,np  

c       write(*,*) spec(i)
c       if (nspec(i).ne.0) then
c         spec(i) = spec(i) ! /float(nspec(i))
c       else
c         spec(i) = 0.               
c       endif

c       write(555,*) i,specu(i),specv(i),specw(i)

c     enddo

      close (555)
      return
      end

C*******************************************************************
      subroutine testftv
C*******************************************************************
C-PD:4/13/03. Tests accuracy of FFT transforms.
C-Performs FT of sample function in x-direction and by inversely FTing
C-checks on accuracy.
      include 'dim.h'

      include 'comwave.h'
      include 'comgridptxyz.h'

      complex rx,expf,sum

      dimension ufin(nx),ufout(nx)
      complex ftfin(nx)
      dimension uf(nx)

      pi2 = 8.0*atan(1.0)
      rx=-1.0*cmplx(0.0,1.0)
C-Set input function values
      do i=1,nx
        xx = x(i)*pi2
        ufin(i) = f0(xx)
      enddo

C-Calculate Fourier coefficients (From Canuto et al. (2.1.22))
      do k=1,nx
        sum = cmplx(0.,0.)
        do j=1,nx
          arg = -float(k)*x(j)*pi2
          expf = cmplx( cos(arg),sin(arg) )
          ampf = cmplx(ufin(j),0)
          sum = sum + ampf*expf          
        enddo
        ftfin(k) = sum/float(nx)
c       write(*,*) k,ftfin(k)
      enddo

C-Recompute initial function
      do j=1,nx
        sum = cmplx(0.,0.)
        do k=1,nx
          arg = float(k)*x(j)*pi2
          expf = cmplx( cos(arg),sin(arg) )
          ampf = ftfin(k)
          sum = sum + ampf*expf          
        enddo
          ufout(j) = real(sum)
      enddo
      
C-Dump out data
c     do j=1,nx
c       write(*,*) j,ufin(j),ufout(j)
c     enddo
 
      return
      end


c
      function f0(x)
c
      pi2 = 8.0*atan(1.0)
      f0 = cos(x)

      return
      end

C**************************************************************
      function ampnoise(ak,akc,akmax,expk,ampmax)
C**************************************************************
C-PD:4/9/03: Noise generation function
C-The amplitude follows a prescribed function k^-11/6 where
C-k is the wavenumber. This corresponds to a total energy
C-of spectrum of k^-5/3 over a spherical surface of radius
C-k in wavenumber space.
C-The phase is selected as a random number between (-pi,pi)

      logical clim

      
      parameter (cnoise=1.e1,tiny = 1.e-30)
	pi2 = 8.0*atan(1.0)
	c1=sqrt(6.)



      clim = ((ak.gt.akc).and.(ak.le.akmax))
      if (clim) then
        ampnoise = ampmax*(ak**expk)/(0.5*pi2*c1)
      else
        ampnoise = 0.
      endif

      return
      end


C****************************************************************
      function umeanprof(umax,rr2)
C****************************************************************
C-7/27/03: Mean profile function for u-velocity
C-umax = maximum defect velocity at centerline.
C-rr2 = radial distance squared
C-r0 = normalization constant
      parameter (r0=0.05)

      umeanprof = umax*(exp(-0.5*rr2/r0**2.))

      return
      end

C****************************************************************
      function rmsprof(umax,rr2)
C****************************************************************
C-PD: 6/29/03: Is function of Gaussian window for rms-profile
C-of noise.
C-umax = maximum defect velocity at centerline.
C-rr2 = radial distance squared
C-r0 = normalization constant

      parameter (a = 0.03375, b=1./0.15, rp = 0.0125, rg = 0.04375)

      rmsprof = umax*a*(b+rr2/rp**2.)*exp(-0.5*rr2/rg**2.)

      return
      end
