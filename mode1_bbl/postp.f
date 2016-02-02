C-Subroutines in this file do various postprocessing of
C-sorts. I.e. energetic, lengthscale calculations.
C-Also contains files associated with relaxation process 
C-(This is intimately related to postprocessing)



C*********************************************************************     
      subroutine postp(initenergy,t,brunt,tstart)    
C*********************************************************************
C-Core subroutine
      include 'dim.h'
 
      include 'comflow.h'
      include 'comenergy.h'
      include 'comenergyw.h'
      include 'comlengthsc.h'
      include 'comsubd.h'

c-Temporarily convert flow back to Physical Space and then
c-back to Fourier-space

      write(*,*) 'POSTPROCESSING DATA'

      call  horfft (u,1,nz)
      call  horfft (v,1,nz)
      call  horfft (w,1,nz)
      call  horfft (temp,1,nz)

C-Evaluate mean u-velocity at all (y,z) by averaging along x-axis
      irelax = 0
      call umeancalc(irelax)
C-Evaluate dumean/dz at all (y,z)
      call dudzmeancalc
C-Evaluate dudymeancalc at all (y,z)
      call dudymeancalc
c-Calculate perturbation values for u at each point in the flow.
      call perturbcalc(irelax)
c-Calculate wake lengthscales (Needs to be done before major
c-postprocessing calculation so that wake center & radii are known).
c     Drop this (due to weird conflicts with compiler flag)
c     call wakelengthcalc
c-Calculate energetic quantities (LHS & RHS of TKE budget)
      call energycalc_postp

c-Convert back to Fourier space
      call  horfft (u,-1,nz)
      call  horfft (v,-1,nz)
      call  horfft (w,-1,nz)
      Call  horfft (temp,-1,nz)

c-Dump data to file
      write(60,'(1x,F10.4,17(2x,E12.6))')
     >(t-tstart),(t-tstart)*brunt,
     >tketot,u1tke,u2tke,u3tke,prody,prodz,prod,
     >bflux,epsilon,epssgs,epstot,
     >chi,chisgs,chitot,
     >tkemean,tempvar

      write(61,'(1x,F10.4,17(2x,E12.6))')
     >(t-tstart),(t-tstart)*brunt,
     >tketotw,u1tkew,u2tkew,u3tkew,prodyw,prodzw,prodw,
     >bfluxw,epsilonw,epssgsw,epstotw,
     >chiw,chisgsw,chitotw,
     >tkemean,tempvarw

      write(62,'(1x,F10.4,8(2x,E12.6))')
     >(t-tstart),(t-tstart)*brunt,
     >(udefy+udefz)/2.,udefy,udefz,wakely,wakelz,yc,zc

      write(*,*) 'Dissipation rate estimate#2',epsilon,epssgs
      return
      end



C*********************************************************************
      subroutine relaxation(initenergy,t,brunt,irelaxend)
C*********************************************************************
C-Core subroutine
      include 'dim.h'

      include 'comflow.h'
      include 'comenergy.h'
      include 'comsubd.h'
c-Temporarily convert flow back to Physical Space and then
c-back to Fourier-space

      write(*,*) 'Performing Relaxation'

      call  horfft (u,1,nz)
      call  horfft (v,1,nz)
      call  horfft (w,1,nz)
      call  horfft (temp,1,nz)

C-Evaluate mean u-velocity at all (y,z) by averaging along x-axis
      irelax = 1 
      call umeancalc(irelax)
C-Evaluate dumean/dz at all (y,z)
      call dudzmeancalc
C-Evaluate dudymeancalc at all (y,z)
      call dudymeancalc

c-Calculate perturbation values for u at each point in the flow.
      call perturbcalc(irelax)

C-Now reset mean and rms perturbation profiles
      call relax_reset(initenergy)

c-Convert back to Fourier space
      call  horfft (u,-1,nz)
      call  horfft (v,-1,nz)
      call  horfft (w,-1,nz)
      call  horfft (temp,-1,nz)

c-Check to see if Production is balanced by
c-sum of dissipation rate and buoyancy flux

c-Dump data to file
C-DATA DUMP NOW DONE IN POSTP
c     write(55,'(1x,F10.4,2x,F10.4,2x,E12.6,2x,E12.6,2x,E12.6,
c    >           2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6)')
c    >t,t*brunt,tketot,prody,prodz,prod,bflux,epsilon,ratio

      return
      end

C*********************************************************************       
      subroutine umeancalc(irelax)
C*********************************************************************
C-Purpose to calculate umean and Tmean across all k=nz x-lines
C-Calculates also v & w mean values (as well as those of SGS tensor)
C-PD-8/2/02: Calculates mean values for v & w used in relaminarization
C-process at top element.
C-irelax: Indicates if this subroutines is called as part of relaxation
C-process.
      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comgridptxyz.h'
      include 'comsgs.h'


C-calculate mean values of u along x direction at every (y,z) line.
C-Initialize     
      do k=1,nz
        do j=1,ny 
          umean(j,k)=0.
          vmean(j,k)=0.
          wmean(j,k)=0.

          t11mean(j,k) = 0.
          t12mean(j,k) = 0.
          t13mean(j,k) = 0.
          t22mean(j,k) = 0.
          t23mean(j,k) = 0.
          t33mean(j,k) = 0.
        enddo
      enddo

      do k=2,nz
         do j=1,ny
            do i=1,nx
               umean(j,k) = umean(j,k) + u(i,j,k)
               vmean(j,k) = vmean(j,k) + v(i,j,k)
               wmean(j,k) = wmean(j,k) + w(i,j,k)

               t11mean(j,k) = t11mean(j,k) + t11(i,j,k)
               t12mean(j,k) = t12mean(j,k) + t12(i,j,k)
               t13mean(j,k) = t13mean(j,k) + t13(i,j,k)
               t22mean(j,k) = t22mean(j,k) + t22(i,j,k)
               t23mean(j,k) = t23mean(j,k) + t23(i,j,k)
               t33mean(j,k) = t33mean(j,k) + t33(i,j,k)           

            enddo
            umean(j,k) = umean(j,k)/nx
            vmean(j,k) = vmean(j,k)/nx
            wmean(j,k) = wmean(j,k)/nx

            t11mean(j,k) = t11mean(j,k)/nx
            t12mean(j,k) = t12mean(j,k)/nx          
            t13mean(j,k) = t13mean(j,k)/nx          
            t22mean(j,k) = t22mean(j,k)/nx          
            t23mean(j,k) = t23mean(j,k)/nx          
            t33mean(j,k) = t33mean(j,k)/nx          

c           if (j.eq.nyh) write(*,*) umean(j,k),vmean(j,k)
          enddo
        enddo

C-Calculate mean value T(z) of temperature. This is used
C-for proper calculation of buoyancy flux.
C-Go to Fourier Space
        call  horfft (temp,-1,nz)

        do k=1,nz
          tempmean(k)=temp(1,1,k)
        enddo
  
        call horfft(temp,1,nz)     


C-This part calculates the Reynolds avgd. part of
C-of the SGS tensor to be used in calculating SGS fluxes
 100    format (1x,f12.6,2x,4(e12.6,2x))

        return
        end




C*********************************************************************
      subroutine dudzmeancalc
C*********************************************************************
C-Purpose to calculate ddz derivative of  umean at all (y,z) positions.
C-This is a modified version of what existed in the global subdomain
C-code. Scans (y,z) profile of velocity in y=const. columns and
C-calculates dudz in each element. Derivative continuity is assumed
C-across element boundaries.
      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comwave.h'
      include 'comleng.h' 

      call dpdzyzplane(umean,dumeandz)

      return
      end



C*********************************************************************
      subroutine dudymeancalc
C*********************************************************************
C-Purpose to calculate ddy derivative of  umean at all (y,z) positions.
      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comfft.h'
   
      dimension sc1((ny+2)*nz)
 
      complex sc1c(ny/2+1,nz)

      equivalence (sc1,sc1c)

C-First copy umean values to a temporary array
      do k=1,nz
        do j=1,ny+2
          k1 = (k-1)*(ny+2) + j
          if (j.le.ny) then
            sc1(k1) = umean(j,k)
          else
            sc1(k1) = 0.
          endif
        enddo
      enddo

C-Convert temporary array to Fourier space
C-U(y,z) is viewed as a 2D array of dimensions nz x (ny+2)
      call fft991 (sc1,work1d,trigsy1d,ifaxy1d,1,ny+2,ny,nz,-1)

C-Calculate derivative in fourier space. Update sc( ) when doing this
C-Careful about defining y-wave# here (Following JAD's advice).
C-This is fully analogous to doing FFT in x-direction in most of other
C-pars of code.
      do k=1,nz
        do j=1,ny/2+1
          k1 = (k-1)*(ny+2) + j      
          ywave = (j - 1)*beta  
          sc1c(j,k) = (0.,1.)*ywave*sc1c(j,k)
        enddo
      enddo

C-Pull out (!! Next time, use prophylactic (stupid Greek humor
C-for anyone who may use this code after me :)- !!) of Fourier space
      call fft991 (sc1,work1d,trigsy1d,ifaxy1d,1,ny+2,ny,nz,1)

C-Copy 1-D array back to 3-D for mean gradient
      do k=1,nz
        do j=1,ny
          k1 = (k-1)*(ny+2) + j
          dumeandy(j,k) = sc1(k1)
        enddo
      enddo

      return
      end



C*********************************************************************
      subroutine perturbcalc(irelax)
C*********************************************************************
C-Calculates perturbation values of u using mean value
C-computed in umean calc. For now u is u'.
C-Do this for v & w in order to reset mean.
C-Not necessary to do this for temperature
C-irelax: Indicates if this subroutines is called as part of relaxation
C-process.
      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comrms.h'

      if (irelax.eq.1) then
        do k=1,nz
          do j=1,ny
            urms(j,k)=0.
            vrms(j,k)=0.
            wrms(j,k)=0.
          enddo
        enddo
        temprms = 0.
      endif

      do k=1,nz
         do j=1,ny
            do i=1,nx
C-Calculate u-velocity perturbations
               u(i,j,k)=u(i,j,k)-umean(j,k)

C-After perturbations have been calculated calculate their rms values
C-(In the initial condition, the only perturbations are with respect
C-to the mean U-velocity). Averaging is performed along homogeneous
C-x-direction (This is only used in the relaxation procedure)
               if (irelax.eq.1) then
                 urms(j,k) = urms(j,k) + u(i,j,k)**2.
                 vrms(j,k) = vrms(j,k) + v(i,j,k)**2.       
                 wrms(j,k) = wrms(j,k) + w(i,j,k)**2.
                 temprms = temprms + temp(i,j,k)**2.
               endif

C-This following line should be used only when relaxation
C-is applied to temperature.
c              temp(ijk)=temp(ijk)-tempmean(k)
             enddo
            
             if (irelax.eq.1) then 
               urms(j,k) = sqrt(urms(j,k)/nx)
               vrms(j,k) = sqrt(vrms(j,k)/nx)
               wrms(j,k) = sqrt(wrms(j,k)/nx)
             endif

           enddo
         enddo

         if (irelax.eq.1) temprms = sqrt(temprms/(nx*ny*nz))

        return
        end
      
      
      
     

C*********************************************************************       
      subroutine umeanreset(zlen)
C*********************************************************************
C-Reset umean values equal to initial values            
C-Set v and w mean values to 0.
      include 'dim.h'

      include 'comwave.h'
      include 'comflow.h'
      include 'commeanv.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'comgridptxyz.h'

      pi2 = 8.0*atan(1.0)
      xl = pi2/alpha
      yl = pi2/beta
C-**PD**: umax is the maximum value of the defect velocity
C-Set through input file now (2/14/02)
c     umax=3.20e-3
c     umax=0.
      z00 = 0.5*zlen
      y00 = 0.5*yl

C-Now reset mean
      write(*,*) 'RESETTING MEAN PROFILE for u-velocity'
      do k=2,nz
         zz=z(k)
         do j=1,ny
            yy= y(j)
            rr2=((zz-z00)**2+(yy-y00)**2)
            u0=umeanprof(umax,rr2)
            umean(j,k) = u0

c           if (j.eq.nyh) write(*,*) yy,zz,umean(j,k)
c           do i=1,nx
C-This following line should be commented out.
C-It's a remnant of trying to avoid temperature solution
C-during relaxation process with time-splitting code.
c              temp(ijk) = tempmean(k)
c           enddo
         enddo
       enddo

 
        return
        end
      

C*****************************************************************
      subroutine energycalc_postp
C*****************************************************************
C-Calculates the following fluctuating TKE budget quantities:
C-
C-tketot =  turbulent kinetic energy.
C-prod = TKE production
C-epsilon = viscous TKE dissipation of resolved scales
C-epssgs = SGS flux of TKE to small scales
C-bflux = Vertical Buoyancy Flux
C-chi = molecular dissipation of temperature variance at resolved scales
C-chisgs = SGS flux of temperature variance to small scales
C-
C-Then outputs data
      logical condwakein

      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'comtkeprev.h'
      include 'comenergy.h'
      include 'comenergyw.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comrms.h'
      include 'comsgs.h'
      include 'comscratch.h'
      include 'comlengthsc.h'
        
C-ALERT: Later insert these arrays into a common block
C-or apply deallocation process.
      include 'dimgrad.h'
      include 'complexflow.h'
      include 'complexgrad.h'

C-Calculate velocity gradients first
      include 'equivflow.h'
      include 'equivgrad.h'

c     dimension tmp(ny,nz)
C-The dthd_i arrays are the spatial gradients of the 
C-fluctuating component of the filtered temperature
C-This is stored in the array ox of the scratch common block. 
C-It has been already calculated through the
C-in subsroutine sgscalc_post. o2 is the array in complex space 
C-that stores this variable.
c     dimension dtdx(nxpp,ny,nz),dtdy(nxpp,ny,nz),dtdz(nxpp,ny,nz)
c     complex dtdxc(nxhp,ny,nz),dtdyc(nxhp,ny,nz),dtdzc(nxhp,ny,nz)
c     equivalence (ox,o1) , (dtdx,dtdxc), (dtdy,dtdyc), 
c    >(dtdz,dtdzc)
      dimension sr(3,3)


C-Convert fluctuating velocities to spectral space
      call  horfft (u,-1,nz)
      call  horfft (v,-1,nz)
      call  horfft (w,-1,nz)
      call  horfft (temp,-1,nz)

      call dpdz(u,dudz)
      call dpdz(v,dvdz)
      call dpdz(w,dwdz)
      call dpdz(temp,dtdz)

c
      do k = 1,nz
       do j = 1,ny
        do i = 1,nxhp
         dudxc(i,j,k) = (0.,1.)*xw(i)*uc(i,j,k)
         dudyc(i,j,k) = (0.,1.)*yw(j)*uc(i,j,k)
         dvdxc(i,j,k) = (0.,1.)*xw(i)*vc(i,j,k)
         dvdyc(i,j,k) = (0.,1.)*yw(j)*vc(i,j,k)
         dwdxc(i,j,k) = (0.,1.)*xw(i)*wc(i,j,k)
         dwdyc(i,j,k) = (0.,1.)*yw(j)*wc(i,j,k)
         dtdxc(i,j,k) = (0.,1.)*xw(i)*tc(i,j,k)
         dtdyc(i,j,k) = (0.,1.)*yw(j)*tc(i,j,k)
        enddo
       enddo
      enddo

C-Convert velocities and gradients back to physical space
      
      call  horfft (u,1,nz)
      call  horfft (v,1,nz)
      call  horfft (w,1,nz)
      call  horfft (temp,1,nz)
      call  horfft (dudx,1,nz)
      call  horfft (dudy,1,nz)
      call  horfft (dudz,1,nz)
      call  horfft (dvdx,1,nz)
      call  horfft (dvdy,1,nz)
      call  horfft (dvdz,1,nz)
      call  horfft (dwdx,1,nz)
      call  horfft (dwdy,1,nz)
      call  horfft (dwdz,1,nz)
      call  horfft (dtdx,1,nz)
      call  horfft (dtdy,1,nz)
      call  horfft (dtdz,1,nz)



C-Calculate TKE budget quantities
C-(Bottom boundary not included)
C-All variables ending in w indicate values calculated
C-strictly within the wake ellipsoid
      tkemean = 0.
      tketot = 0.
      u1tke = 0.
      u2tke = 0.
      u3tke = 0.
      prody = 0.
      prodz = 0.
      prod = 0.
      epsilon = 0.
      epssgs = 0.
      bflux = 0.
      tempvar = 0.
      chi = 0. 
      chisgs = 0.

      tketotw = 0.
      u1tkew = 0.
      u2tkew = 0.
      u3tkew = 0.
      prodyw = 0.
      prodzw = 0.
      prodw = 0.
      epsilonw = 0.
      epssgsw = 0.
      bfluxw = 0.
      tempvarw = 0.
      chiw = 0.
      chisgsw = 0.
     
      nyc = nyh+1
      nzc = (nz-1)/2+1 
      
      do k=2,nz
         do j=1,ny
            do i=1,nx

C-Is this point inside the wake ellipsoid ?
               yy2 = ((y(j)-yc) / wakely)**2.
               zz2 = ((z(k)-zc) / wakelz)**2.
               condwakein = ((zz2.le.1.).and.(yy2.le.1.))

c              if (yy2.le.1) write(*,*)
c    >         z(k),zc,wakelz,zz2
c              if ((zz2.le.1).and.(yy2.le.1.)) write(*,*)
c              if (condwakein) write(*,*)
c    >         y(j),yc,wakely,yy2

c              write(*,*) condwakein

               tketot = tketot +
     +         (u(i,j,k)**2.+v(i,j,k)**2.+w(i,j,k)**2.)
 
               u1tke=  u1tke + u(i,j,k)**2.
               u2tke = u2tke + v(i,j,k)**2.
               u3tke = u3tke + w(i,j,k)**2.

C-Production has components in both y & z directions:
               prody = prody  
     >               -dumeandy(j,k)*(u(i,j,k)*v(i,j,k))
               prodz = prodz
     >               -dumeandz(j,k)*(u(i,j,k)*w(i,j,k))

c              temp(ijk) = 0.
               bflux = bflux + w(i,j,k)*temp(i,j,k)           
              
               tempvar = tempvar + temp(i,j,k)**2. 
C-Viscous dissipation of TKE at resolved scales
               epsilon = epsilon +
     >         (dudx(i,j,k)**2.)+(dudy(i,j,k)**2.)+(dudz(i,j,k)**2.) +
     >         (dvdx(i,j,k)**2.)+(dvdy(i,j,k)**2.)+(dvdz(i,j,k)**2.) +      
     >         (dwdx(i,j,k)**2.)+(dwdy(i,j,k)**2.)+(dwdz(i,j,k)**2.)

C-Calculate the above quantities for a region only inside the wake
C-ellipsoid
               if (condwakein) then

                 tketotw = tketotw +
     +           (u(i,j,k)**2.+v(i,j,k)**2.+w(i,j,k)**2.)

                 u1tkew = u1tkew + u(i,j,k)**2.
                 u2tkew = u2tkew + v(i,j,k)**2.
                 u3tkew = u3tkew + w(i,j,k)**2.

C-Production has components in both y & z directions:
                 prodyw = prodyw
     >                 -dumeandy(j,k)*(u(i,j,k)*v(i,j,k))
                 prodzw = prodzw
     >                 -dumeandz(j,k)*(u(i,j,k)*w(i,j,k))

                 bfluxw = bfluxw + w(i,j,k)*temp(i,j,k)

                 tempvarw = tempvarw+temp(i,j,k)**2.

C-Viscous dissipation of TKE at resolved scales
                 epsilonw = epsilonw +
     >           (dudx(i,j,k)**2.)+(dudy(i,j,k)**2.)+(dudz(i,j,k)**2.) +
     >           (dvdx(i,j,k)**2.)+(dvdy(i,j,k)**2.)+(dvdz(i,j,k)**2.) +
     >           (dwdx(i,j,k)**2.)+(dwdy(i,j,k)**2.)+(dwdz(i,j,k)**2.)
 
               endif
                 
 

C-SGS TKE flux to unresolved scales 
C-First calculate:
C-a) strain rate tensor for unresolved scales 
               sr(1,1) = dudx(i,j,k)
               sr(2,2) = dvdy(i,j,k)
               sr(3,3) = dwdz(i,j,k)
               sr(1,2) = 0.5*(dudy(i,j,k)+dvdx(i,j,k))
               sr(1,3) = 0.5*(dudz(i,j,k)+dwdx(i,j,k))
               sr(2,3) = 0.5*(dvdz(i,j,k)+dwdy(i,j,k))
C-b) perturbations to Reynolds avgd. value of SGS tensor
               t11fluct = t11(i,j,k) - t11mean(j,k)
               t12fluct = t12(i,j,k) - t12mean(j,k)
               t13fluct = t13(i,j,k) - t13mean(j,k)
               t22fluct = t22(i,j,k) - t22mean(j,k)
               t23fluct = t23(i,j,k) - t23mean(j,k)
               t33fluct = t33(i,j,k) - t33mean(j,k)

                
               sgs =    t11fluct*sr(1,1) +
     >                  t22fluct*sr(2,2) +
     >                  t33fluct*sr(3,3) +
     >               2.*t12fluct*sr(1,2) +
     >               2.*t13fluct*sr(1,3) +
     >               2.*t23fluct*sr(2,3)

               epssgs = epssgs - sgs

C-Specifically for wake interior
               if (condwakein) then

                 epssgsw = epssgsw - sgs 

               endif

C-Molecular dissipation of temp. var. at resovled scales
               chi = chi +
     >               (dtdx(i,j,k)**2.)+
     >               (dtdy(i,j,k)**2.)+
     >               (dtdz(i,j,k)**2.)

C-SGS temp. var.  flux to unresolved scales
               sgs = th1(i,j,k)*dtdx(i,j,k) +
     >               th2(i,j,k)*dtdy(i,j,k) +
     >               th3(i,j,k)*dtdz(i,j,k)

               chisgs = chisgs - sgs

C-Calculate these last two quantities if point falls inside
C-wake ellipsoid

              if (condwakein) then
               
                chiw = chiw +
     >                 (dtdx(i,j,k)**2.)+
     >                 (dtdy(i,j,k)**2.)+
     >                 (dtdz(i,j,k)**2.)

C-SGS temp. var.  flux to unresolved scales
                 sgs = th1(i,j,k)*dtdx(i,j,k) +
     >                 th2(i,j,k)*dtdy(i,j,k) +
     >                 th3(i,j,k)*dtdz(i,j,k)

                 chisgsw = chisgsw - sgs 
                   
              endif

            enddo


c-Calculate mean TKE
            tkemean = tkemean + umean(j,k)**2.
          enddo
        enddo

        tketot = tketot*0.5
        tkemean = tkemean*0.5*float(nx)
        bflux = bflux*ra*pr 
c       write(*,*) bflux,ra
        epsilon = epsilon*pr
        epstot = epsilon + epssgs
        chi = chi*xkappa
        chitot = chisgs + chi
        prod = prody+prodz

        tketotw = tketotw*0.5
        bfluxw = bfluxw*ra*pr
c       write(*,*) bflux,ra
        epsilonw = epsilonw*pr
        epstotw = epsilonw + epssgsw
        chiw = chiw*xkappa
        chitotw = chisgsw + chiw
        prodw = prodyw+prodzw

c-Add back to u the contribution from the mean
        do k=2,nz
          do j=1,ny
            do i=1,nx
               u(i,j,k) = u(i,j,k)+umean(j,k)
            enddo
          enddo
        enddo

        write(*,*) 'DISSIPATION RATES',epssgs,epsilon
        write(*,*) 'Scalar DRs',chisgs,chi
c       write(*,*) chitotw,chiw,chisgsw
        return
        end



C*********************************************************************
        subroutine relax_reset(initenergy)
C*********************************************************************
C-Resets mean velocity profile as well as that of x-averaged rms
C-value of perturbation velocities to mantain both constant profile
C-but also total rms value of perturbation field.
        include 'dim.h'

        parameter(irelaxtemp = 0,temprmspre=0.125)
        logical condtemp

        include 'comflow.h'
        include 'commeanv.h'
        include 'comgridptxyz.h'
        include 'comsubd.h'
        include 'comrms.h'
        
        

        call umeanreset(zlen)

C************************************************************
C-This is from old procedure where the constant which multiplied
C-TKE values to keep them constant was globally the same
C-New procedure discussed in updated version of Dommermuth et al.
C-requires local coefficient as f unction of (y,z)
C************************************************************
C-Compare with TKE level at previous timestep
C-and update accordingly. Also add mean velocity to u-field.
C
C**r = ratio of TKE at time n over that at time n-1
C**np = number of points for which TKE needs to be updated

C
C-Update velocity at each point by
C sqrt(2/(np*r))
c       np = (nz-1)*ny*nx
C         if (initenergy.eq.0) then
C           r = tketot/tkeprev
C         write(*,*) 'r,tketot,tkeprev,initenergy'
C         write(*,*) r,tketot,tkeprev,initenergy
C         else
C           r = 1.
C         endif

C         factke = sqrt(1./r)

          utke = 0.
          vtke = 0.
          wtke = 0.
          temprmsnew = 0.

c         do j=1,ny
c           do k=1,nz
c             tmp(j,k) = 0.
c           enddo
c         enddo

C-Check and see if rms temperature fluctuations exceed
C-prescribed level
          condtemp = ((irelaxtemp.eq.1).and.(temprms.ge.temprmspre)) 
          write(*,*) 'Temperature Relaxation: ',condtemp
          if (condtemp) facrmst = temprmspre/temprms 
           
          do k=2,nz
           do j=1,ny

C-Find multiplicative factors for u,v & w that will ensure rms of velocity
C-at this location will remain constant
            if (initenergy.eq.0) then
              facrmsu = urms0(j,k)/urms(j,k)
              facrmsv = vrms0(j,k)/vrms(j,k)
              facrmsw = wrms0(j,k)/wrms(j,k)
            else
              facrmsu = 1.
              facrmsv = 1.
              facrmsw = 1.
            endif

            do i=1,nx

                 u(i,j,k) = u(i,j,k)*facrmsu
                 v(i,j,k) = v(i,j,k)*facrmsv
                 w(i,j,k) = w(i,j,k)*facrmsw

c                tempfluct = temp(i,j,k) - tempmean(k)

C-Modulate temperature fluctuation if it exceeds limits
C-Then add it to mean profile (We assume this doesn't get mixed during relaxation,
C-i.e. it stays constant)
                 if (condtemp) then
c                  tempfluct = temp(i,j,k) - tempmean(k)
                   temp(i,j,k) = temp(i,j,k)*facrmst
                   temprmsnew = temprmsnew + temp(i,j,k)**2.
                   temp(i,j,k) = tempmean(k) + temp(i,j,k)
                 endif
                  
                 utke = utke + u(i,j,k)**2.
                 vtke = vtke + v(i,j,k)**2.
                 wtke = wtke + w(i,j,k)**2.
 
                 tketot = utke + vtke + wtke

C-Add back mean velocity profile which has already been reset
                 u(i,j,k) = u(i,j,k)+umean(j,k)

c                tmp(j,k) = tmp(j,k) + w(i,j,k)**2.            
            enddo
            
            urms0(j,k) = urms(j,k)*facrmsu
            vrms0(j,k) = vrms(j,k)*facrmsv
            wrms0(j,k) = wrms(j,k)*facrmsw

c           tmp(j,k) = sqrt(tmp(j,k)/nx)

           enddo
          enddo

          if (condtemp) temprmsnew = sqrt(temprmsnew/(nx*ny*nz))


c         write(*,*) 'TKE components'
c         write(*,*) utke,vtke,wtke
c         write(*,*) 'Old t. rms, new t. rms, factor'
c         write(*,*) temprms,temprmsnew,facrmst

c         j = nyh+1
c         do k=1,nz
c           write(*,*) z(k),tmp(j,k)
c         enddo

          tketot = tketot*0.5
C-Set initenergy = 0 for subsequent steps
          initenergy = 0
C
C--------------------------------------------------------------------------
C--------------------------------------------------------------------------

        return
        end



c*********************************************************************
      subroutine wakelengthcalc
c*********************************************************************
C-Calculates the spanwise and vertical width of the mean wake profile
C-Normally this routine should use some least-squares Gaussian fit.
C-This uses the refined version obtained from Geoff Spedding in August 02.
C
C-Step 1: Estimate fit quantities by performing manual calculation.
C-Wake thickness is assumed to be at location 0.1 of defect velocity centerline
C-value.
C
C-Step 2: Apply curve-fit to obtain more exact estimate
C-For now we simply move from the centerline a distance such that
C-the mean velocity magnitude is 0.1 its centerline value.
      include 'dim.h'

C-Parameter nn must be compatible with the value used in subroutine curvefit.
      parameter (facwidth=0.1,nn=1000)
     
      include 'comwave.h'
      include 'comlengthsc.h' 
      include 'commeanv.h'
      include 'comleng.h'
      include 'comenergy.h'
      include 'comsubd.h'
      include 'comgridptxyz.h'

      dimension utmpy(nn),utmpz(nn),ctmpy(nn),ctmpz(nn)

C-Calculate maximum defect velocity 
C-This is a calculation designed to produce a first estimate for
C-the maximum defect velocity and the wake width to immediately
C-insert into the curvefit routine.
      nyc = nyh+1
      nzc = (nz-1)/2+1
      udef0 = umean(nyc,nzc)

C-Focus on centerplane (Oxy)
      k=nzc
C
C-Scan the mean profile in y-direction to calculate L_y
      do j=1,nyc-1
        if ((umean(j,k)/udef0).ge.facwidth) then 
          jmin = j
          goto 10
        endif
      enddo

  10  continue

      do j=ny,nyc,-1 
        if ((umean(j,k)/udef0).ge.facwidth) then 
          jmax = j
          goto 20
        endif
      enddo

   20 wakely0 = 0.5*abs(y(jmax)-y(jmin))

C-Focus on centerplane (Oxz)
      j=nyc
C-Scan the mean profile in z-direction to calculate L_z
      do k=1,nzc-1     
        if ((umean(j,k)/udef0).ge.facwidth) then
          kmin = k
          goto 30
        endif
      enddo

  30  continue

      do k=nz,nzc,-1
        if ((umean(j,k)/udef0).ge.facwidth) then
          kmax = k
          goto 40
        endif
      enddo
 

   40 wakelz0 = 0.5*abs(z(kmax)-z(kmin))
 
     
C-Now use curve-fit to get more accurate estimate of fits
      do j=1,ny
        utmpy(j) = umean(j,nzc)
        ctmpy(j) = y(j)
      enddo
      call curvefit(udef0,y(nyc),wakely0,utmpy,
     >ctmpy,ny,udefy,yc,wakely)


      do k=1,nz
        utmpz(k) = umean(nyc,k)
        ctmpz(k) = z(k)
      enddo
      call curvefit(udef0,z(nzc),wakelz0,utmpz,
     >ctmpz,nz,udefz,zc,wakelz)

c     write(*,*) '---'
c     write(*,*) udefy,yc,wakely
c     write(*,*) udefz,zc,wakelz
      return
      end

c*********************************************************************
      subroutine dpdzyzplane (a,b)
c*********************************************************************
C-Modification of dpdz. Calculates dpdz for data available only
C-On a yz-plane. This is a different version than that which is
C-found in old global subdomain code.
      include 'dim.h'
      include 'comsubd.h'
      include 'comleng.h'
c     common /index/  li(ntot),lj(ntot),lk(ntot)
c
      dimension a(ny,nz),b(ny,nz)
      dimension tmp1r(nzloc),dtmp1r(nzloc)
c

C-Scan over all y-values
      do 100 j=1,ny
C-Scan over each subdomain in z-direction
        do 200 ks=1,nsubd
          
          ybar = zh(ks)
C-Transfer data to temporary array
          do  kloc=1,nzloc
            k = (ks-1)*nzloc + kloc
            tmp1r(kloc) = a(j,k)   
          enddo

C-Now calculate derivative
          call dpdzcol(tmp1r,dtmp1r,d,zpts,nzm,ybar,1,nzloc)


C-Transfer data back to final array
          do  kloc=1,nzloc
            k = (ks-1)*nzloc + kloc
            b(j,k) = dtmp1r(kloc)   
          enddo        

 200    continue          
 100  continue

      return
      end      

