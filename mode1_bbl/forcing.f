C******************************************************
C-This file contains all forcing functions associated
C-with soliton at top free surface
C******************************************************

C-CORRECTION (PD: 2/27/04): Wave propagates from left
C-to right with phase velocity umax (C in notes).
C-A -umax factor should appear in expressions for 
C-density and pressure for the wave.

!
!-PD-Lyon-6/30/2015: Modified forcing functions to
!-represent MODE-1 wave
!


C******************************************************
C-Soliton wavenumber
C******************************************************

      function solwvn(zeta0,zlen,rh1h2)
      real h1,h2 

      h2 = zlen
      h1 = rh1h2*h2

      coeff = 0.75*abs(h1-h2)/(h1*h1*h2*h2)
     
      solwvn = sqrt(zeta0*coeff)
      return
      end



C*****************************************************
C-Boundary Layer Thickness (Scaling estimate)
C*****************************************************  

      function deltabl(umax,zlen,xnu,zeta0,rh1h2,solks)

      reinv = 1./sqrt(umax*zlen/xnu)
      f1 = (1.-rh1h2)/(zeta0/zlen)
      f2 = sqrt(0.5/(solks*zlen)) 
      f3 = sqrt(1./rh1h2)

      deltabl = reinv*f1*f2*f3*zlen
      

c     deltabl = sqrt(xnu*solks*zlen/(umax*zeta0))/(0.5/solks)

      return
      end



C*****************************************************
C-Horizontal dependence of soliton, function A(x) (see
C-Larry's notes).
C-Derivatives and integrals
C*****************************************************

      function asol(x,xlen,solks,a0)
      xs = x -1./2.*xlen
      sech2 = 1./(cosh(solks*xs)**2.)

      asol = a0*sech2

      
      return
      end

cc

      function dasoldx(x,xlen,solks,a0)
      xs = x -1./2.*xlen
      sech2 = 1./(cosh(solks*xs)**2.)

      dasoldx =  -2.*a0*solks*sech2*tanh(solks*xs)

      return
      end

cc

      function d2asoldx2(x,xlen,solks,a0)
      xs = x -1./2.*xlen
      sech2 = 1./(cosh(solks*xs)**2.)

      d2asoldx2 =  2.*a0*(solks**2.)*sech2*
     >             (2. * tanh(solks*xs)**2. -
     >              sech2)
     

      return
      end

cc

      function d3asoldx3(x,xlen,solks,a0)
      xs = x -1./2.*xlen
      sech2 = 1./(cosh(solks*xs)**2.)

      d3asoldx3 =  8.*a0*(solks**3.)*
     >             sech2*tanh(solks*xs)*
     >             ( 2.*sech2 - tanh(solks*xs)**2. )

      return
      end      

C***********************************************************
C-Mean Density Gradient


C***********************************************************
       function rhograd(z,zlen,grav,rho0,brunt)

        rhograd = -rho0*(bruntz(z,zlen,brunt)**2.)/grav
 
       return
       end
      
       function rhomean_unif(z,zlen,grav,rho0,brunt)
C-Assume for all practical purposes that the density is 1
C-at the bottom of the domain.

      rhomean_unif = 1.+ rhograd(z,zlen,grav,rho0,brunt)*z

      return
      end

C***********************************************************
C-Mean Density Profile (ARBITRARY STRAT.)
C***********************************************************
        function rhomean(z,zlen,grav,rho0,brunt)
C-Assume for all practical purposes that the density is 1
C-at the bottom of the domain.

      rhomean = 1.+ rhograd(z,zlen,grav,rho0,brunt)*z

      return
      end


      function bruntz_unif(z,brunt)

      bruntz_unif = brunt

      return
      end



C***********************************************************
C-Brunt-Vaisala Frequency Vertical Profile (ARBITRARY STRAT.)
C***********************************************************
      function bruntz(z,zlen,brunt)

C-Dt represents non-dimensional thickness of thermocline
       bruntz = brunt

      return
      end

C***********************************************************
C-2nd Vertical derivative of Mean Density Profile (ARBITRARY STRAT.)
C***********************************************************
      function d2rhomeandz2(z,zlen,grav,rho0,brunt)
!  C-Dt represents non-dimensional thickness of thermocline
!      parameter (dt=0.15, a=0.98611, b=0.01389, dn=1.)
!      z0 = 0.2625

!      arg = 4.*(z-z0)/(dt*0.45)
!      arg0 = 4./(dt*0.45)

!      sech2 = 1./(cosh(arg)**2.)
!      tanh2 = 1. - dn*tanh(arg)

!      d2rhomeandz2 = -brunt**2.*rho0/grav*arg0*sech2 *
!     >               (-2.*a*tanh(arg) - b*dn)
      d2rhomeandz2 = 0.0 

      return
      end

C***********************************************************
C-Non-Dimensional Velocity Shear Profile
C***********************************************************
C*** This part is to be modified when we use the Gaussian shear
C*** profile..A.M.A  03-04-2008
      function fshear(z,zdim,zsh,xacu)
       implicit none
       real z,zdim,xacu,fshear,zsh,dt,z0,zdom,gam
       parameter (dt = 0.13, z0 = 0.2625, zdom = 0.45)

       gam = 4.0/(dt*zdom)
       fshear = 0.5*(tanh(gam*(z - z0))+1) 

       return
      end


C***********************************************************
C-Non-Dimensional Velocity Shear Profile
C***********************************************************
      function dfsheardz(z,zdim,zsh,xacu)
      implicit none
      real dfsheardz,zsh,xacu,z,zdim,fshear,dt,z0,zdom,gam
      parameter (dt = 0.13, z0 = 0.2625, zdom = 0.45)

      gam = 4.0/(dt*zdom)
      dfsheardz = 0.5*gam*(1.0/cosh(gam*(z - z0)))**2

      return
      end


C-INFORMATION ASSOCIATED WITH "VIRTUAL PADDLE" OF
C-SLINN AND RILEY (1998).
C-This part of the code contains forcing functions
C-for the u & w velocity fields and the density perturbation.

C- Second derivative of the shear flow A.M.A. 03-04-2008

      function d2fsheardz(z,zdim,zsh,xacu)
      implicit none
      real d2fsheardz,xacu,z,zdim,zsh,dfsheardz,fshear,dt,z0,zdom,gam,
     >     hyptan,hypsec
      parameter (dt = 0.13, z0 = 0.2625, zdom = 0.45)

      gam = 4.0/(dt*zdom)
      hyptan = tanh(gam*(z - z0))
      hypsec = 1.0/cosh(gam*(z - z0))
      d2fsheardz = -1.0*(gam**2)*hyptan*hypsec**2

      return
      end


C*************************************************************
c-Localization function and its first derivative
C*************************************************************
      function fzloc(z,zcen,zdim,xbcrit)
         implicit none 
         real fzloc,xbcrit,z,zcen,zdim,PI
        PI=4.0*atan(1.0)
         fzloc = exp(-xbcrit*(z-zcen)*(z-zcen)/(zdim**2))
c         if ((z <= zcen+0.5*PI).and.(z >= zcen-0.5*PI)) then
c         fzloc=(sin((z-zcen-0.5*PI)/zdim))**2
c         else
c           fzloc=0.0
c         endif 
      return
      end

           function fxgaus(x,xcen,zdim,xfrac,p)
         implicit none
         real fxgaus,sigx,x,xcen,zdim,xfrac
         integer p,PI
             sigx=xfrac*zdim
             fxgaus=exp(-0.5*((x-xcen)/sigx)**p)
        return
        end

      function fzgaus(z,zcen,zdim,xbcrit,p)
         implicit none
         real fzgaus,xbcrit,z,zcen,zdim,PI,Zlo,Zth
         integer p
 
          PI=4.0*atan(1.0)
c          Zth=2.0*zdim*sqrt(-log(0.1)/xbcrit)
c          Zlo=zcen-(zdim*sqrt(-log(0.1)/xbcrit))
          fzgaus=exp(-xbcrit*((z-zcen)/zdim)**p)
c sine function 
c          if((z.gt.Zlo).and.(z.lt.(Zlo+Zth)))then
c           fzgaus=sin((z-Zlo)*PI/Zth)
c          else
c           fzgaus=0.0   
c           endif
        return
        end

c

      function fzlocp(z,zcen,zdim,xbcrit)
        implicit none
        real fzlocp,fzloc,xbcrit,z,zcen,zdim,PI
         PI=4.0*atan(1.0)
     
        fzlocp = -2.*xbcrit*(z-zcen)*fzloc(z,zcen,zdim,xbcrit)
     >            /(zdim**2)
c         if ((z <= zcen+0.5*PI).and.(z >= zcen-0.5*PI)) then
c         fzlocp=(2./zdim)*(sin((z-zcen-0.5*PI)/zdim))*
c     >           (cos((z-zcen-0.5*PI)/zdim))
c         else
c         fzlocp=0.0
c         endif 
      return
      end

C*************************************************************
c-Forcing functions which generate IW field
C*************************************************************
! 
!
!-PD-Lyon-7/1/15: I cleaned up these functions to correspond
!-to a MODE-1 wave and eliminate some of Ammar's chaos.
!
!
!
c-Remember: these are the appropriately dimensionalized versions
c-of what is found in eqs. (2.29) thru (2.32) of Slinn & Riley
c
!-----------
c-u-momentum
!-----------

c zdim is the vertical length scale used for dimensionalizing
c the nondimensional source terms in Slinn & Riley A.M.A 02-14-2008  
c xkcrit,xmcrit , and xomega  are dimensional quantities 
!
!-PD-2/7/15: All forcing functions need to be divided through with
!-zdim/umax .
!-The magnitude of the eigenfunction is given by xacrit*umax
!
      function forceu(x,z,t,
     >                xacrit,xkcrit,xmcrit,xomega,
     >                zdim,umax)
c- rampu-up of source terms A.M.A. 03-03-2008 
           implicit none 
           real sinf,cosf,veldim,forceund,wzp
           real x,z,t,brunt,PI,phit,forceu,
     >                xacrit,xkcrit,xmcrit,xomega,
     >                zdim,umax

!-Target velocity is u(x,z,t) = (-1/k)*W'(z)*sin(k*x-omega*t)

       wzp = xmcrit*cos(xmcrit*z) ! Vertical Eigenfunction z-derivative
       sinf = sin( xkcrit*x - xomega*t )
       cosf = cos( xkcrit*x - xomega*t ) 

       forceu = -(1/xkcrit)*
     >          (xacrit*umax)*  
     >          wzp*sinf*
     >          umax/zdim*
     >          phit(t,xomega)


        return
        end

!-----------
c-w-momentum
!-----------

      function forcew(x,z,t,
     >                xacrit,xkcrit,xmcrit,xomega,
     >                zdim,umax)
         implicit none
          real sinf,cosf,veldim,forcewnd,wz
          real x,z,t,brunt,phit,forcew,
     >                xacrit,xkcrit,xmcrit,xomega,
     >                zdim,umax

!-Target velocity is w(x,z,t) = W(z)*cos(k*x-omega*t)

        wz = sin(xmcrit*z) ! Vertical Eigenfunction
        cosf = cos( xkcrit*x  - xomega*t )

        forcew = (xacrit*umax)*
     >           wz*cosf*
     >           (umax/zdim)*
     >           phit(t,xomega)

      return
      end
       
!---------------------
c-Density perturbation
!---------------------

      function forcerho(x,z,t,
     >                  xacrit,xkcrit,xmcrit,xomega,
     >                  zdim,umax,
     >                  brunt,rho0,grav)

         implicit none 
         real sinf,cosf,veldim,forcerhond,wz
         real xomegandim, rhodim,forcerho,rho0
         real x,z,t,fzlocp,brunt,grav,
     >             phit, xacrit,xkcrit,xmcrit,xomega,
     >              PI,zdim,umax

!-Target density is rho(x,z,t) = -(1/omega)*N^2*(rho0/g)*W(z)*sin(k*x-omega*t)

        PI=4.0*atan(1.0)
        sinf = sin( xkcrit*x  - xomega*t )

        wz = sin(xmcrit*z) ! Vertical Eigenfunction
        
        forcerho = (xacrit*umax)*
     >             wz*sinf*
     >             -(1/xomega)*brunt*brunt*(rho0/grav)*
     >            (umax/zdim)*
     >            phit(t,xomega)

        return
        end

!----------------------------------------
!-- Temporal envelope of forcing function
!----------------------------------------


         function phit(t,xomega)
         implicit none
         real t,xomega,phit,twave,pi2,omegaenv

         pi2 = 8.*atan(1.0)
         twave = pi2/xomega
         omegaenv = xomega/4.

!-Forcing is turned off after one wave period

         if (t <= 8.*twave) then
           phit = 1. !cos(omegaenv*t)
         else
           phit = 0.
         endif

         return
         end

         function indexfi(zfi)
           include 'dim.h'
           include 'comgridptxyz.h'
        
                       do k=1,nz
               if(z(k).ge.zfi) then
                indexfi=k
                exit
               endif
               enddo
           return
           end
!!!!!!!!!!!!!!! HORIZONTAL AVERAGE !!!!!!!!!!!!!!1
         function horav(Ufield,kk,p)
          include'dim.h'
          REAL, DIMENSION(nxpp,ny,nz) :: Ufield
           integer :: kk,p
           real ::sumtemp,horav
            sumtemp=0.0
            do i=1,nx
             sumtemp=sumtemp+(Ufield(i,1,kk)**p)
            enddo
           horav=sumtemp/nx
         return
          end
       function horav2(Ufiel,Wfiel,kk)
       !!this is so damn expensive on 512 grid because I define two huge arrays!
          include'dim.h'
          REAL, DIMENSION(nxpp,ny,nz) :: Ufield,Wfield
           integer :: kk
           real ::sumtemp,horav2
            sumtemp=0.0
            do i=1,nx
             sumtemp=sumtemp+(Ufiel(i,1,kk)*Wfiel(i,1,kk))
            enddo
           horav2=sumtemp/nx
         return
          end

!!!!!!!!!!!!!!!!!!! VERTICAL AVERAGE !!!!!!!!!!!!!!!!!!!!!!
          function verav(Vfield)
!--------------------------------------------------
           include'dim.h'
           include 'comsubd.h'
!--------------------------------------------------
          REAL, DIMENSION(nz) :: Vfield
           real, dimension(1:nzloc) ::aloc
           real :: ybar ,outloc,sumztemp,verav
           integer :: ks,kloc,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              sumztemp=0.0
           do ks=1,nsubd   
              ybar = zh(ks)
               aloc=0.0
               outloc=0.0
        do kloc = 1,nzloc
           k = (ks-1)*nzloc + kloc
           aloc(kloc)=Vfield(k)
        enddo
        call integrcol(aloc,outloc,ndim,ybar)
            sumztemp = sumztemp + outloc
         enddo
           verav=sumztemp/zlen

          return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

          function phiz(x,z,xkcrit,xmcrit,zcen,zdim,xbcrit,xlen)
          
c         implicit none
C***  Le,Ue are the lower and upper limits of the z Rayleigh absorber
C** Phi(z)/A =0 at z=Le, =1 at z=Ue and dph/dz=0 at z=Le
c         common /domlength/  ZL,Zcur  !  domain height
c         parameter(rzo=.25)! 20% of Zlen/2 is dedciated to the absorber layer
c         parameter(zAmp=1.0) !Amplitude of phi(z)
c         real z,ao,a1,a2,Cz,Le,Ue,zAmp,rzo,ZL,phiz,Dab,zp,Zcur
c         logical FLAGZ
cc******Defining the vertical limits of the upper and bottom absorbing layers******
c        if (z >0.5*ZL) then
c        Le=(1.0-0.5*rzo)*ZL
c        Ue=ZL
c        FLAGZ=.TRUE.
c        else
c        Le=0.5*rzo*ZL ! this is the upper boundary of the bottom absorber layer
c        Ue=0.0        ! the lower boundary of the bottom abosrber layer 
c        FLAGZ=.FALSE.
c        endif
c**********Evaluating polynomial coefficients ****************       
         
c              Cz=Le**2-2.0*Le*Ue+Ue**2
c              a2=1.0/Cz  
c              a1=-2*a2*Le
c              ao=a2*Le**2                    
c          if(FLAGZ.and.z.ge.Le) phiz=Amp*(ao+a1*z+a2*z**2)
c          if(FLAGZ.and.z.lt.Le) phiz=0.0               
c          if(.NOT.FLAGZ.and.z.le.Le) phiz=Amp*(ao+a1*z+a2*z**2)
c          if(.NOT.FLAGZ.and.z.gt.Le) phiz=0.0
c translate coordinates to a new system in which the two absorbers have
c the same coordinate z'=|z-zo| where zo is half the domain height
c          Dab=abs(0.5*ZL*(1-rzo))
c          zp=abs(z-0.5*ZL)
c          if(zp.lt.Dab) then
c
c           zp=(z-0.5*ZL)
c           if(zp.lt.Dab) then
c          phiz=0.0
c          else
c          Le=Dab
c          Ue=abs(0.5*ZL)
c          Cz=Le**2-2.0*Le*Ue+Ue**2
c          a2=1.0/Cz
c          a1=-2*a2*Le
c          ao=a2*Le**2
c          phiz=zAmp*(ao+a1*zp+a2*zp**2)          
c          endif
c
         

c         return
c         end

           implicit none

           real phiz,x,z,xkcrit,xmcrit,ampx,ampz,spot,
     >          rdxwall,lspongedx,rspongedx,spongedz,
     >          xlen,pi2,lcorner,rcorner,zcen,zdim,xbcrit,
     >          locamp,vreduce,hreduce,l,
     >          fxgaus,fzloc,reduce,onoff,phizlim

           parameter (ampx= 1. ,ampz= 1., locamp= 3.0, onoff = 1.)

      if(onoff.eq.1.) then

!        if(z.gt.0.175) then 

!          phiz = 0.0

!        else

           pi2= 8.0*atan(1.0)

           rspongedx= 1.2*abs(pi2/xkcrit)
           lspongedx= 0.4*rspongedx
           spongedz= 0.7*abs(pi2/xmcrit)

           rdxwall= xlen - x
           
           lcorner= (spongedz/lspongedx)*x - z
           rcorner= (spongedz/rspongedx)*rdxwall - z

           if(rdxwall.lt.rspongedx) then
                if(rcorner.lt.0) then
!                     l= sqrt((2*(z - 0.1225))**2 + rdxwall**2)
!                   if((l.lt.rspongedx).and.(z.gt.0.1225))then
!                     spot= (l - rspongedx)/rspongedx
!                     phiz= ampx*(sin(0.25*pi2*spot))**2
!                   elseif(z.lt.0.1225)then
                     spot= (rdxwall-rspongedx)/rspongedx
                     phiz= ampx*(sin(0.25*pi2*spot))**2
!                   else
!                     phiz= 0.0
!                   endif
                else
                     spot= (z-spongedz)/spongedz
                     phiz= ampz*(sin(0.25*pi2*spot))**2
                endif
                
c           vreduce= fzloc(z,zcen,zdim,xbcrit)
c           hreduce= fxgaus(x,1.5,2.64,1./4.2919,2)
c           reduce= locamp*(vreduce * hreduce)
c 
c                if(reduce.lt.phiz) then
c                phiz= phiz-reduce
c                elseif(reduce.gt.phiz) then
c                phiz= 0.0
c                endif
                
           elseif(x.lt.lspongedx) then
                if(lcorner.lt.0) then
!                     l= sqrt((z-0.1225)**2 + x**2)
!                   if((l.lt.lspongedx).and.(z.gt.0.1225)) then
!                     spot= (l - lspongedx)/lspongedx
!                     phiz= ampx*(sin(0.25*pi2*spot))**2
!                   elseif(z.lt.0.1225)then
                     spot= (x-lspongedx)/lspongedx
                     phiz= ampx*(sin(0.25*pi2*spot))**2
!                   else
!                     phiz = 0.0
!                   endif
                else
                     spot= (z-spongedz)/spongedz
                     phiz= ampz*(sin(0.25*pi2*spot))**2
                endif

c           vreduce= fzloc(z,zcen,zdim,xbcrit)
c           hreduce= fxgaus(x,1.5,2.64,1./4.2919,2)
c           reduce= locamp*(vreduce * hreduce)
c
c                if(reduce.lt.phiz) then
c                phiz= phiz-reduce
c                elseif(reduce.gt.phiz) then
c                phiz= 0.0
c                endif
 
           elseif(z.lt.spongedz) then
                spot= (z-spongedz)/spongedz
                phiz= ampz*(sin(0.25*pi2*spot))**2

c           vreduce= fzloc(z,zcen,zdim,xbcrit)
c           hreduce= fxgaus(x,1.5,2.64,1./4.2919,2)
c           reduce= locamp*(vreduce * hreduce)
c 
c                if(reduce.lt.phiz) then
c                phiz= phiz-reduce
c                elseif(reduce.gt.phiz) then
c                phiz= 0.0
c                endif
            
           else
                phiz= 0.0

!           endif

         endif

       endif

      end
