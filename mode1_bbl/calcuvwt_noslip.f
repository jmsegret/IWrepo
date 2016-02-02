C*************************************************************************
      subroutine calcuvwt(izmode,ispecfilter,irelax,
     >cpenspecfilter)
C*************************************************************************

c-Calculates All velocity Components
c-Remember ! We are starting out in spectral space.
C-**PD: 2/15/03**: Derivatives calculated directly
C-with Costa and Don formulae.
C-**PD: 1/3/03**: This routine consists of the viscous/
C-pressure step treatment and any necessary updates.
C-I.e. fractional stepping is enhanced with spectral
C-multidomain formulation.
      parameter(tiny = 1.e-16,ifilterp=1,ifilterw=1)
      include 'parambctype.h'
      include 'paramvarfilterleg.h'
      include 'dim.h'
 
      include 'comflow.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comcoeff.h'
      include 'commatrix.h'
      include 'comzeta.h'
      include 'comwfunc.h'
      include 'comddzwfunc.h'
      include 'combc.h'
      include 'comtime.h'
      include 'comparamphys.h'
      include 'comsubd.h'
      include 'comhomogstore.h'
      include 'comgridptxyz.h'
      include 'comrotpressbc.h'
      include 'comtimeadvance.h'
      include 'comcritlayer.h'


      complex rx

      logical lzeromode,cpenspecfilter,modecheck
c-All these are auxilliary  functions. I'd rather have more
c-1-D functions like this than store the result into the same
c-array over and over again !

      dimension rhsinr(nz),rhsini(nz)
      dimension fvr(nzloc),fvi(nzloc),foutr(nzloc),fouti(nzloc)
      dimension fr(nz),fi(nz) !,dfr(nz),dfi(nz),d2fr(nz),d2fi(nz) 

      dimension ub(nxpp,ny,nz),vb(nxpp,ny,nz),wb(nxpp,ny,nz)
      complex ucol(nz),vcol(nz),wcol(nz),divvel(nz)
      complex dudx(nz),dvdy(nz),dwdz(nz)
c     complex dudxold(nz),dvdyold(nz),dwdzold(nz)
c     complex dudxnew(nz),dvdynew(nz),dwdznew(nz)
c     dimension divrold(nz),diviold(nz),divrnew(nz),divinew(nz)
      complex phi(nxhp,ny,nz)

c     dimension rhsinrlc(nzloc),rhsinilc(nzloc)

c     dimension solprtr(nsubd,nzloc),solprti(nsubd,nzloc),
c    >          solhom1r(nsubd,nzloc),
c    >          solhom2r(nsubd,nzloc)


        common /Misc/tzmode


      parameter (niter=500)
         real  PI
         PI=4.0*atan(1.0)

      rx = -1.0*cmplx(0.,1.)

      open(250,file='testin.dat')
c     write(250,*) 'VARIABLES = "z",
c    >              "Re(dudx)-old","Im(dudx)-old",
c    >              "Re(dvdy)-old","Im(dvdy)-old",
c    >              "Re(dwdz)-old","Im(dwdz)-old",
c    >              "Re(Div)-old","Im(Div)-old",
c    >              "Re(dudx)-new","Im(dudx)-new",
c    >              "Re(dvdy)-new","Im(dvdy)-new", 
c    >              "Re(dwdz)-new","Im(dwdz)-new",
c    >              "Re(Div)-new","Im(Div)-new"'
     
c

c-Scan over all wavenumbers
      do 100 j=1,ny 
         ii = 0
         do 200 i=1,nxpp,2 
           ii = ii + 1

           lzeromode = ((ii.eq.1).and.(j.eq.1)) 
           modecheck = ((ii.eq.1).and.(j.eq.1))
           ak = xsq(ii) + ysq(j)

C***************
C*PRESSURE STEP*
C***************

C-Now solve equation for Surrogate Pressure Phi.
C-See KIO: This is the averaged pressure over (n)->(n+1)
C-interval.

C-Careful with zero mode
           if (lzeromode) then

             do k=1,nz
               phi(ii,j,k) = cmplx(0.,0.)
             enddo

           else

C-Set up beta coefficient
             coeffbeta = ak
             epsilon = 1./coeffbeta
             call penalty_setup(epsilon)

C-Nulify divvel(...) array and calculate grad(u^*)
C-(i.e. divergence of intermediate velocity field)

             do k=1,nz
                  divvel(k) = cmplx(0.,0.)
                  ucol(k) = cmplx(su(i,j,k),su(i+1,j,k))
                  vcol(k) = cmplx(sv(i,j,k),sv(i+1,j,k))
                  wcol(k) = cmplx(sw(i,j,k),sw(i+1,j,k))
             enddo

c
c-Store in temporary variables calculate divergence
c-of velocity for this column of data.
             xwin =  xw(ii)
             ywin =  yw(j)
             call divcolc(ucol,vcol,wcol,xwin,ywin,divvel,
     >                    dudx,dvdy,dwdz)

C-Set up RHS (imaginary and real)
             do k=1,nz


               rhsinr(k) = real(divvel(k))/dt
               rhsini(k) = real(rx*divvel(k))/dt

             enddo

C-Set up BC's. Work with homogeneous BC's or KIO
C-recommended improved BC

             ifuncreal = 0.
             ifbctop = ineumann
c            ifbctop = idirichlet
             ifbcbot = ineumann
c            ifbcbot = idirichlet
             bcrtop =   wnl(i,j,2) + wviscrot(i,j,2)
             bcrbot =   wnl(i,j,1) + wviscrot(i,j,1)
             bcitop =   wnl(i+1,j,2) + wviscrot(i+1,j,2)
             bcibot =   wnl(i+1,j,1) + wviscrot(i+1,j,1)

             if (modecheck) write(*,*) 'Pressure Bcs',
     >       bcrtop,bcitop,bcrbot,bcibot,wnl(i,j,2),wnl(i+1,j,2)

C-Now solve for surrogate pressure PHI

C-Scan over all subdomains and calculate
C-solutions following penalty method
C-(Penalty method works fine provided epsilon is finite)

c           write(*,*) 'Examining Pressure for Mode',ii,j
C-Call penalty method for calculation of solution
C-(Method borrowed from Hesthaven, 97)
             call penalty_main(epsilon,ifuncreal,
     >                       bcrbot,bcibot,bcrtop,bcitop,
     >                       ifbcbot,ifbctop,rhsinr,rhsini,
     >                       fr,fi)

C-Filter spectrally legendre modes in each subdomain
             if (cpenspecfilter)  then
c            if (ifilterp.eq.1) then

               do ks=1,nsubd

                 ybar = zh(ks)
C-Update temporary input array
                 do kloc=1,nzloc
                   k = (ks-1)*nzloc + kloc
                   fvr(kloc) = fr(k)
                   fvi(kloc) = fi(k)
                 enddo
C-Filter in Legendre spectral space
                 call filter_legenpen(fvr,foutr,ybar,iprt,ks)
                 call filter_legenpen(fvi,fouti,ybar,iprt,ks)
C-Transfer data from temporary output array
                 do kloc=1,nzloc
                   k = (ks-1)*nzloc + kloc
c                if (lzeromode) write(*,*)
c    >           z(k),fvr(kloc),foutr(kloc)
                   fr(k) = foutr(kloc)
                   fi(k) = fouti(kloc)
                 enddo

               enddo

             endif

C-PD:5/5/03(Brown): Follow Jan's suggestion. Average at interfaces
C-Immediately after solution of Helmholtz eqn.
c            call interface_avg(fr)
c            call interface_avg(fi)

C-Now update PHI array
             do k=1,nz
               phi(ii,j,k) = cmplx(fr(k),fi(k))
c              if (modecheck) write(250,'(1x,5(e12.6,2x))')
c    >z(k),
c    >    -rhsinr(k)/epsilon,-rhsini(k)/epsilon,
c    >fr(k),fi(k)
             enddo
c-Finished with non-zero mode

          endif

 200    continue
 100  continue
 
c     pause
c
C-Now apply PRESSURE CORRECTION to the velocity
C-ub,vb,wb is the output (incompressible velocity)
c
      call update(ub,vb,wb,u,v,w,phi,dt)

C-NOTE: NO NEED TO CALL EXBC. No intermediate BC's used
C-except for pressure BC

c     goto 500


C************************
C-V I S C O U S  S T E P*      
C************************
c       write(*,*)'xbcrit',xbcrit,'zcen',zcen
   
cxbcrit=15.72,
c    zcen=5.5,
c    xacu=108.3,
c    zsh=3.85 $

      do 125 j=1,ny 
        ii = 0
        do 225 i=1,nxpp,2 
          ii = ii + 1

           lzeromode = ((ii.eq.1).and.(j.eq.1))
           modecheck = ((ii.eq.2).and.(j.eq.1))
           ak = xsq(ii) + ysq(j)

C*************
C*U-VELOCITY**
C*************
C*************** Nullifying the zero mode for t < 2*PI/xomega ****
C *************** A.M.A. May,28th, 2008 ***************
c               write(*,*)'KLi',KLi,KUi
               
c
c               if(t.lt.tzmode)izmode=1
c               if(t.ge.tzmode)izmode=0
c                if (lzeromode) then !.and.t < izmode*6000.0*PI/xomega) then
c                       do k=KLi,KUi
c                       fr(k) = 0.
c                       fi(k) = 0.
cc                       end do
c                
c
c                else 
 
C-Set up beta coefficient
           coeffbeta = ak + bdfv*facvel
           epsilon = 1./coeffbeta
           call penalty_setup(epsilon)
        
C-Set up RHS (imaginary and real)
           do k=1,nz

             rhsinr(k) = -ub(i,j,k)/(xnu*dt)
             rhsini(k) = -ub(i+1,j,k)/(xnu*dt)

           enddo


C-Set up BC's. This is not a real function
C-dudz = 0 at top, u = 0 at bottom
           ifuncreal = 0.
           ifbctop = ineumann
!          ifbctop = idirichlet
!          ifbcbot = ineumann
           ifbcbot = idirichlet
           bcrtop = 0.! uw(i,j,2)
           bcrbot = 0. ! uwuw(i,j,1)
           bcitop = 0. !uw(i+1,j,2)
           bcibot = 0. !uw(i+1,j,1)

C-Call penalty method for calculation of solution
C-(Method borrowed from Hesthaven, 97)
           call penalty_main(epsilon,ifuncreal,
     >                       bcrbot,bcibot,bcrtop,bcitop,
     >                       ifbcbot,ifbctop,rhsinr,rhsini,
     >                       fr,fi)

C-Filter spectrally legendre modes in each subdomain
           if (cpenspecfilter)  then

             do ks=1,nsubd

               ybar = zh(ks)
C-Update temporary input array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
                 fvr(kloc) = fr(k)
                 fvi(kloc) = fi(k)
               enddo
C-Filter in Legendre spectral space
               call filter_legenpen(fvr,foutr,ybar,iprt,ks)
               call filter_legenpen(fvi,fouti,ybar,iprt,ks)
C-Transfer data from temporary output array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
                 fr(k) = foutr(kloc)
                 fi(k) = fouti(kloc)

               enddo

             enddo

           endif

C***************END OF NULLIFYING ZERO MODE LOOP , A.M.A.05-28-2008**
c           endif 
!check z_grid for z(k) in the forcing region
*fzloc(z,zcen,zdim,xbcrit)
c               write(*,*)z(nz),zcen,zdim,xbcrit,
c     >  fzloc(z(nz),zcen,zdim,xbcrit)
C************************************************************
c              if (lzeromode) then
               if (lzeromode) then !.and.(t.gt.2.0*PI/xomega)) then   
C*********************************************************************
!!! Note the width parameter xbcrit is real however the exponent p is integer
                     do k=1,nz
c             fr(k) =(1.0-fzgaus(z(k),zcen,zdim,9.6,2))*fr(k)
             fr(k) =(1.0-fzgaus(z(k),zcen,zdim,xbcrit,2))*fr(k)
c     >                  (1.0-exp(-0.5*(t/100.0)**2))
             fi(k) =(1.0-fzgaus(z(k),zcen,zdim,xbcrit,2))*fi(k)
c     >                  (1.0-exp(-0.5*(t/100.0)**2))
c             write(*,*)1-fzgaus(z(k),zcen,zdim,16.,2),z(k)
                         enddo
cccc fzgaus(z,zcen,zdim,xbcrit,p)
c                            do  k=1,nz    !KLi,KUi
c            fr(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fr(k)
c            fi(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fi(k)
c                           end do
c                     endif
C****************************************************************
               endif
C***************************************************************                  
 
C-Now calculate actual intermediate velocity u^*
           do k=1,nz 
             u(i,j,k) = fr(k) 
             u(i+1,j,k) = fi(k) 
           enddo

C*************
C*V-VELOCITY**
C*************
C*************** Nullifying the zero mode for t < 2*PI/xomega ****
C *************** A.M.A. May,28th, 2008 ***************
c                    write(*,*)'Tzmode=',tzmode
c               if(t.lt.tzmode)izmode=1
c               if(t.ge.tzmode)izmode=0
c
c                if (lzeromode) then!.and.t < izmode*6000.0*PI/xomega) then
c                       do k=KLi,KUi
c                         fr(k) = 0.
c                         fi(k) = 0.
c                       end do
c                else

C-Set up beta coefficient
           coeffbeta = ak + bdfv*facvel
           epsilon = 1./coeffbeta
           call penalty_setup(epsilon)

C-Set up RHS (imaginary and real)
           do k=1,nz
        
             rhsinr(k) = -vb(i,j,k)/(xnu*dt)
             rhsini(k) = -vb(i+1,j,k)/(xnu*dt)
   
           enddo

C-Set up BC's. This is not a real function
C-dvdz = 0 at top, v = 0 at bottom
           ifuncreal = 0.
           ifbctop = ineumann
c          ifbctop = idirichlet
           ifbcbot = idirichlet
!          ifbcbot = ineumann
           bcrtop = 0.         
           bcrbot = 0.         
           bcitop = 0.         
           bcibot = 0.                    

C-Call penalty method for calculation of solution
C-(Method borrowed from Hesthaven, 97)
           call penalty_main(epsilon,ifuncreal,
     >                       bcrbot,bcibot,bcrtop,bcitop,
     >                       ifbcbot,ifbctop,rhsinr,rhsini,
     >                       fr,fi)

C-Filter spectrally legendre modes in each subdomain
           if (cpenspecfilter)  then
 
             if (modecheck) iprt = 1
             do ks=1,nsubd

               ybar = zh(ks)
C-Update temporary input array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
                 fvr(kloc) = fr(k)
                 fvi(kloc) = fi(k)
               enddo
C-Filter in Legendre spectral space
               call filter_legenpen(fvr,foutr,ybar,iprt,ks)
               call filter_legenpen(fvi,fouti,ybar,iprt,ks)
C-Transfer data from temporary output array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
c                if (lzeromode) write(*,*)
c    >           z(k),fvr(kloc),foutr(kloc)
                 fr(k) = foutr(kloc)
                 fi(k) = fouti(kloc)
               enddo

             enddo

           endif

C***************END OF NULLIFYING ZERO MODE LOOP , A.M.A.05-28-2008**             
c           endif
C****************************************
c                     if (lzeromode) then
c               if (lzeromode) then !.and.(t.gt.2.0*PI/xomega)) then   
c                     do k=1,nz
c             fr(k) =(1.0-fzgaus(z(k),zcen,zdim,9.6,2))*fr(k)
cc     >                  (1.0-exp(-0.5*(t/100.0)**2))
c             fi(k) =(1.0-fzgaus(z(k),zcen,zdim,9.6,2))*fi(k)
cc     >                  (1.0-exp(-0.5*(t/100.0)**2))
cc              write(*,*)(1.0-fzgaus(z(k),zcen,zdim,108.3,2))
c               write(*,*)(1.0-exp(-0.5*(t/100.0)**2))
c                         enddo

C*********************************************************************
c                      if(t.lt.8.0*PI/xomega) then
c                        do k=1,nz    !KLi,KUi
c              fr(k) = 0.0     !(1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fr(k)
c              fi(k) = 0.0 !(1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fi(k)
c                        end do
c c                     else
c                            do  k=1,nz    !KLi,KUi
ccc            fr(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fr(k)
c            fi(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fi(k)
cic                           end do
c                     endif
C****************************************************************
c               endif 

C*******************************************

C-Now calculate actual intermediate velocity u^*
           do k=1,nz
             v(i,j,k) = fr(k) 
             v(i+1,j,k) = fi(k) 
           enddo

C************
C*W-VELOCITY*
C************

!          goto 1500
C-CAREFUL: If this is the zero Fourier mode (0,0)
C-set up corresponding w-velocity to zero
           if (lzeromode) then

             do k=1,nz
               fr(k) = 0.
               fi(k) = 0.
c              wb(1,1,k)=0.
c              wb(2,1,k)=0.
             end do

C-Treat all other Fourier modes now


           else


C-Set up beta coefficient
             coeffbeta = ak + bdfv*facvel
             epsilon = 1./coeffbeta
             call penalty_setup(epsilon)

C-Set up RHS (imaginary and real)
             do k=1,nz

               rhsinr(k) = -wb(i,j,k)/(xnu*dt)   
               rhsini(k) = -wb(i+1,j,k)/(xnu*dt) 
!              if (modecheck) write(*,*) rhsinr(k),rhsini(k)
 
             enddo


C-Set up BC's. This is not a real function
C-dudz = 0 at top, u = 0 at bottom
             ifuncreal = 0.
c            ifbctop = ineumann
           ifbctop = idirichlet
c          ifbcbot = ineumann
             ifbcbot = idirichlet
             bcrtop = 0.         
             bcrbot = 0.         
             bcitop = 0.         
             bcibot = 0.         

C-Call penalty method for calculation of solution
C-(Method borrowed from Hesthaven, 97)
             call penalty_main(epsilon,ifuncreal,
     >                       bcrbot,bcibot,bcrtop,bcitop,
     >                       ifbcbot,ifbctop,rhsinr,rhsini,
     >                       fr,fi)
 
C-Filter spectrally legendre modes in each subdomain
             if (cpenspecfilter)  then
c            if (ifilterw.eq.1) then

               do ks=1,nsubd

                 ybar = zh(ks)
C-Update temporary input array
                 do kloc=1,nzloc
                   k = (ks-1)*nzloc + kloc
                   fvr(kloc) = fr(k)
                   fvi(kloc) = fi(k)
                 enddo
C-Filter in Legendre spectral space
                 call filter_legenpen(fvr,foutr,ybar,iprt,ks)
                 call filter_legenpen(fvi,fouti,ybar,iprt,ks)
C-Transfer data from temporary output array
                 do kloc=1,nzloc
                   k = (ks-1)*nzloc + kloc
!                if (modecheck) write(*,*)
!    >           z(k),fvr(kloc),foutr(kloc)
                   fr(k) = foutr(kloc)
                   fi(k) = fouti(kloc)

c                if ((ii.eq.1).and.(j.eq.1)) write(250,*)
c    >       z(k),-rhsinr(k)/epsilon,fr(k)
                 enddo 
               enddo 

             endif

           endif
c      open(653,file='testcalcvar') 
c      write(653,*)i,j,t,xomega,PI,2*PI/xomega


C-Now calculate actual intermediate velocity u^*
           do k=1,nz
             w(i,j,k) = fr(k) 
             w(i+1,j,k) = fi(k) 
           enddo

 1500  continue
c          if (ii.eq.1) write(*,*) ii,j,wb(i,j,1),wb(i+1,j,1)

 225    continue
 125  continue

C*************
C*TEMPERATURE*
C*************

C-SKIP TEMPERATURE CALCULATION WHEN PERFORMING RELAXATION
C-PROCESS (With new relaxation process which accounts for development
C-of temperature field, drop this)
C
C-
C-

C-Rescan over all wavenumbers
C

c     goto 500

      do 150 j=1,ny 
        ii = 0
        do 250 i=1,nxpp,2 
          ii = ii + 1
!        write(*,*) 'Solving Temperature for wavenumbers',i,j

           lzeromode = ((ii.eq.1).and.(j.eq.1))
           modecheck = ((ii.eq.2).and.(j.eq.2))
           ak = xsq(ii) + ysq(j)

C-In this case we are solving for temperature at level n+1/2.
C*************** Nullifying the zero mode for t < 2*PI/xomega ****
C *************** A.M.A. May,28th, 2008 ***************
c               if(t.lt.tzmode)izmode=1
c               if(t.ge.tzmode)izmode=0
c              write(*,*)'KLi',KLi,KUi

c                if (lzeromode) then !.and.t < izmode*6000.0*PI/xomega) then
c                         do k=KLi,KUi
c                          fr(k) = 0.
c                          fi(k) = 0.
c                          end do
c                else

C-Set up beta coefficient
           coeffbeta = ak + bdfv*fact
           epsilon = 1./coeffbeta
           call penalty_setup(epsilon)

C-Set up RHS (imaginary and real)
           do k=1,nz

             rhsinr(k) = - temp(i,j,k)/(dt*xkappa) 
             rhsini(k) = - temp(i+1,j,k)/(dt*xkappa) 

           enddo

C-Set up BC's. This is not a real function
C-Fixed temperatures at both boundaries (Set in Init.f).
           ifuncreal = 0.
!          ifbctop = idirichlet
!          ifbcbot = idirichlet
           
           ifbctop = ineumann
           ifbcbot = ineumann
 
           bcrtop = 0. !tubc(i,j)
           bcrbot = 0. !tbbc(i,j)
           bcitop = 0. !tubc(i+1,j)
           bcibot = 0. !tbbc(i+1,j) 


C-Call penalty method for calculation of solution
C-(Method borrowed from Hesthaven, 97)
           call penalty_main(epsilon,ifuncreal,
     >                       bcrbot,bcibot,bcrtop,bcitop,
     >                       ifbcbot,ifbctop,rhsinr,rhsini,
     >                       fr,fi)

C-Filter spectrally legendre modes in each subdomain
           if (cpenspecfilter)  then

             do ks=1,nsubd

               ybar = zh(ks)
C-Update temporary input array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
                 fvr(kloc) = fr(k)
                 fvi(kloc) = fi(k)
               enddo
C-Filter in Legendre spectral space
               call filter_legenpen(fvr,foutr,ybar,iprt,ks)
               call filter_legenpen(fvi,fouti,ybar,iprt,ks)
C-Transfer data from temporary output array
               do kloc=1,nzloc
                 k = (ks-1)*nzloc + kloc
                 fr(k) = foutr(kloc)
                 fi(k) = fouti(kloc)
               enddo

             enddo

           endif

C-Update temperature array at time level n+1. Careful since result
C-of helmholtz solver is at time n+1/2
C***************END OF NULLIFYING ZERO MODE LOOP , A.M.A.05-28-2008**
c            endif
C**************************************************

c              if (lzeromode) then
c              if (lzeromode) then !.and.(t.gt.2.0*PI/xomega)) then   
cC*********************************************************************
c                     do k=1,nz
c             fr(k) =(1.0-fzgaus(z(k),zcen,zdim,15.8,2))*fr(k)
cc     >                  (1.0-exp(-0.5*(t/100.0)**2))
c c            fi(k) =(1.0-fzgaus(z(k),zcen,zdim,15.8,2))*fi(k)
c     >                  (1.0-exp(-0.5*(t/100.0)**2))
c                         enddo
c                      if(t.lt.8.0*PI/xomega) then
c                        do k=1,nz    !KLi,KUi
c              fr(k) = 0.0     !(1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fr(k)
c              fi(k) = 0.0 !(1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fi(k)
c                        end do
c                      else
c                            do  k=1,nz    !KLi,KUi
cc            fr(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fr(k)
c            fi(k) = (1.0-fzgaus(z(k),zcen,zdim,108.3,2))*fi(k)
c                           end do
c                     endif
C****************************************************************
c               endif
c            fr(k) =(1.0-fzgaus(z(k),zcen,zdim,15.72,2))*fr(k)
c     >                  (1.0-exp(-0.5*(t/100.0)**2))
c            fi(k) =(1.0-fzgaus(z(k),zcen,zdim,15.72,2))*fi(k)
c     >                  (1.0-exp(-0.5*(t/100.0)**2))
c             enddo
c             endif
          

C****************************************************

           do k=1,nz
             temp(i,j,k) = fr(k) 
             temp(i+1,j,k) = fi(k)
           enddo


 500    continue

 250    continue
 150  continue 

c     close(250)
             
      return
      end



c*********************************************************************
      subroutine exbc(phi,uw,vw,ww)
c*********************************************************************
c
c   extrapolate the velocity boundary condition
c   assuming the wall velocity is zero for the lower boundary (1)
c   and the upper boundary (nz) is stress free
c   **PD-03/03/03**: Modify this to handle more general boundary conditions.
c   PD: IMPROVEMENT DUE VERY SOON: Dphi/Dz is only needed at 
c   bottom and top subdomains
c
      include 'dim.h'

      include 'comleng.h'
      include 'comtime.h'
      include 'comwave.h'
      include 'comscratch.h'

c     complex uw(nxhp,ny,*),vw(nxhp,ny,*),ww(nxhp,ny,*),unit,
c    1        o2(nxhp,ny,nz),o1(nxhp,ny,nz)
c     dimension phi(nxpp,ny,nz)
c     equivalence (ox,o1), (oy,o2)
C-CAREFUL: In contrast to update, u,v & w are treated as real
C-variables here.
      complex uw(nxhp,ny,2),vw(nxhp,ny,2),ww(nxhp,ny,2)
      complex phi(nxhp,ny,nz),dphidz(nxhp,ny,nz)
      complex unit,tmp,rx

      rx = -1.0*cmplx(0.,1.)
c
c-Calculate the globall the vertical derivative
      call dpdz (phi,dphidz)

c
      unit = (0.,1.)
      do j=1,ny
       ii = 0
       do i=1,nxpp,2
       ii = ii+1

C-U-velocity: bottom boundary
        tmp = phi(ii,j,1)
c       if (ii.eq.1) write(*,*) j,phi(ii,j,1)
        uw(ii,j,1) = unit*xw(ii)*tmp*dt
C-U-velocity: top boundary
        tmp = dphidz(ii,j,nz)
        uw(ii,j,2) = unit*xw(ii)*tmp*dt

c       if (ii.eq.1)
c    >write(*,*) ii,j,xw(ii),phi(ii,j,1),dphidz(ii,j,nz)
c      >tmp,uw(ii,j,1),uw(ii,j,2)
C-V-velocity: bottom boundary
        tmp = phi(ii,j,1)
        vw(ii,j,1) = unit*yw(j)*tmp*dt
C-V-velocity: top boundary
        tmp = dphidz(ii,j,nz)
        vw(ii,j,2) = unit*yw(j)*tmp*dt

C-W-velocity: bottom boundary
        tmp = dphidz(ii,j,1)
        ww(ii,j,1) = tmp*dt
C-W-velocity: top boundary
        tmp = dphidz(ii,j,nz)
        ww(ii,j,2) = tmp*dt
     
       enddo
      enddo
c
       return
       end



 
c
c*********************************************************************
      subroutine update(u,v,w,ub,vb,wb,phi,dt)
c*********************************************************************
c
c   update the velocity field
c   (**PD: 3/3/03**): May need to make this more flexible
c   to accommodate for a variety of BC's.
c   Could also be modified to avoid introducing all these
c   new complex arrays.
c
      include 'dim.h'
c
      include 'comleng.h'
      include 'comscratch.h'
      include 'comwave.h'
      include 'comgridptxyz.h'
      logical lzeromode,modecheck

      complex phi(nxhp,ny,nz),dphidz(nxhp,ny,nz)

C-NOTE: u,v,w are treated as complex variables here
      complex u(nxhp,ny,nz),v(nxhp,ny,nz),w(nxhp,ny,nz),
     1        ub(nxhp,ny,nz),vb(nxhp,ny,nz),wb(nxhp,ny,nz)
      complex tmp,rx

      rx = -1.0*cmplx(0.,1.)
c
c-Calculate the global vertical derivative
      call dpdz (phi,dphidz)
c
c       do k=2,nz-1  ! for both rigid surface!j
       do k=2,nz    ! for upper surface (nz) free
!      do k=1,nz ! Both surfaces are free slip (PD-Lyon-7/1/15)
        do j=1,ny
         ii = 0
         do i=1,nxpp,2
         ii = ii+1

           modecheck = ((ii.eq.1).and.(j.eq.1))
           lzeromode = ((ii.eq.1).and.(j.eq.1))

c          if (modecheck) write(250,*)
c          if (modecheck) write(*,*)
c    >     z(k),phi(i,j,k),dphidz(i,j,k)
C-Update U and V velocities first
           tmp = phi(ii,j,k)
           u(ii,j,k) = ub(ii,j,k)  - cmplx(dt,0.)*cmplx(0.,xw(ii))*tmp
           v(ii,j,k) = vb(ii,j,k)  - cmplx(dt,0.)*cmplx(0.,yw(j))*tmp

C-Now update W velocity
C-Note that for top surface, w is left unchanged.
C-(Should be zero right ? Answer this question ...)
C-Careful with zero mode: It's always left equal to zero
C-Note: (12/6/03: W should be updated at top free surface too).
           if (lzeromode) then
             w(ii,j,k) = cmplx(0.,0.)
           else
             if (k.lt.nz+1) then
               tmp = dphidz(ii,j,k)
               w(ii,j,k) = wb(ii,j,k) - cmplx(dt,0.)*tmp
             else
               w(ii,j,k) =  wb(ii,j,k)
             endif
           endif

c          if (modecheck) write(*,*) z(k),w(ii,j,k),wb(ii,j,k)
c          if (modecheck) write(250,'(1x,5(e12.6,2x))')
c    >z(k),dummy,dummy,real(u(ii,j,k)),real(rx*u(ii,j,k))
c    >z(k),dummy,dummy,real(v(ii,j,k)),real(rx*v(ii,j,k))
c    >z(k),dummy,dummy,real(w(ii,j,k)),real(rx*w(ii,j,k))

         enddo
        enddo
       enddo

       return
       end           

           
c
c*********************************************************************
      subroutine update_decoy(u,v,w,ub,vb,wb,phi,dt)
c*********************************************************************
c
c   update the velocity field
c   (**PD: 3/3/03**): May need to make this more flexible
c   to accommodate for a variety of BC's.
c   Could also be modified to avoid introducing all these
c   new complex arrays.
c
      include 'dim.h'
c
      include 'comleng.h'
      include 'comscratch.h'
      include 'comwave.h'
      include 'comgridptxyz.h'
      logical modecheck

      dimension phi(nxpp,ny,nz),dphidz(nxpp,ny,nz)

C-NOTE: u,v,w are treated as complex variables here
      complex u(nxhp,ny,nz),v(nxhp,ny,nz),w(nxhp,ny,nz),
     1        ub(nxhp,ny,nz),vb(nxhp,ny,nz),wb(nxhp,ny,nz)
      complex tmp
c
c-Calculate the global vertical derivative
      call dpdz (phi,dphidz)
c
c       do k=2,nz-1  ! for both rigid surfaces
       do k=2,nz    ! for upper surface (nz) free
        do j=1,ny
         ii = 0
         do i=1,nxpp,2
         ii = ii+1

           modecheck = ((ii.eq.24).and.(j.eq.48))

c          if (modecheck) write(250,*)
c          if (modecheck) write(*,*)
c    >     z(k),phi(i,j,k),dphidz(i,j,k)
C-Update U and V velocities first
           tmp = cmplx( phi(i,j,k),phi(i+1,j,k) )
           u(ii,j,k) = ub(ii,j,k) - 0. ! cmplx(0.,1.)*dt*xw(ii)*tmp
           v(ii,j,k) = vb(ii,j,k) - 0. ! cmplx(0.,1.)*dt*yw(j)*tmp

C-Now update W velocity
C-Note that for top surface, w is left unchanged.
C-(Should be zero right ? Answer this question ...)
           if (k.lt.nz) then
             tmp = cmplx( dphidz(i,j,k),dphidz(i+1,j,k) )
             w(ii,j,k) = wb(ii,j,k) - 0. ! cmplx(dt,0.)*tmp
           else
             w(ii,j,k) = w(ii,j,k)
           endif

         enddo
        enddo
       enddo

       return
       end           


c*****************************************************
        subroutine addVelocityNoise(istep,xlen)
c*****************************************************

C Adds noise to the velocity field to induce Kelvin-Helmholz 
C instability. Designed to work in the presence of a fshear function
c which can be found in forcing.f


      include 'dim.h'
      include 'comflow.h'
      include 'comgridptxyz.h'
        real a, b, k0, h0, unoise, wnoise, z0, u0,pi,xlen
        integer nTimes, coin, istep

c h0,a,b and k0 are based on smythe paper

        h0=0.03375
        z0=0.2625

c k0 ALWAYS equals this number (0.8976) for 1st mode disturbances.
        k0= 0.8976

        a=0.05
        b=0.4177
        x0=xlen/2.
        
c u0 is the total change in ushear across the vertical domain

        u0 = 0.015

c horfft gives us the velocity fields in physical space

      call  horfft (u,1,nz)
      call  horfft (w,1,nz)

c This opens 'pert.dat' to write data for a contour plot of our perterbtion field

        open(2020,file='pert.dat')
        write(2020,*) 'ZONE F=POINT, I=', nx, 'J=', nz

c do for everywhere x and z

       do k=1,nz
c        do j=1,ny
          do i=1,nx

c These perterbation functions are also from the smythe paper.
c Unoise and Vnoise must agree with continuity del(u_field) = 0
c 	unoise and vnoise are the noise added.

c		Calculates vertical velocity noise...
                 argx = (x(i)-x0)/h0
                 argz = (z(k)-z0)/h0
                 tnh = tanh(2.*argz)
                 sch = 1./cosh(2.*argz)
                 cs1 = cos(2.*k0*argx)
                 cs2 = cos(k0*argx)
                 unoise = (0.5*u0/k0)*(-cs1 + 2.*b*cs2)*tnh*sch

c		Calculates horizontal velocity noise...
                 sn1 = sin(2.*k0*argx)
                 sn2 = sin(k0*argx)
                 wnoise = 0.5*u0*(sn1 - b*sn2)*sch
               
c		Stores the perterbations in pert.dat...
                 write(2020,*) x(i),z(k),unoise,wnoise

c Noise terms are added to velocity field in physical space.
c They are about a*u0 in magnitude.
                 u(i,1,k)= u(i,1,k)+a*unoise
                 w(i,1,k)= w(i,1,k)+a*wnoise

          enddo
c        enddo
       enddo

       close(2020)

c This takes the new velocity fields with perterbations added and
c transforms them back to fourier space.

      call  horfft (u,-1,nz)
      call  horfft (w,-1,nz)

      end


c*******************************************************************
c              Paul's random number generator July 29, 2010
c       Not used for anything!
c       Varies from -1 to 1.
c       Randomizes with respect to x,y,z and time. Hence the i, j, k,
c 	and istep inputs. 

      function rndm(i,j,k,istep)
         
         implicit none
         real rndm, pi, zz, modi, modj, modk, arg,ri,rj,rk,ristep
         integer i,j,k,istep

C Take the real part of our interger inputs so we combine them with other
C real numbers and do modulo functions. i.e. use amod instead of mod.

         ri = real(i)
         rj = real(j)
         rk = real(k)
         ristep = real(istep)

c Arg is makes a really big number that changes rapidly with z, since
c z is small and in the denominator. Rndm takes the modulo of it, and we don't really
c know where rndm will end up.

         pi = 4.0*atan(1.)
         zz = k*0.0035
         modi = amod(ri,31.)
         modj = amod(rj,17.)
         modk = amod((rk + pi*modi),23.)
         arg=10000.*(amod(ristep,10.)*modi/modj*modk + 5.)/(zz+pi)
         rndm = (2./(2.*pi))*(amod(arg,(2.*pi))) - 1.

c Sometimes the amod function does not work correctly in rndm, and we get numbers
c with a magnitude greater than one. If the number is bigger than one, it goes through
c the loop below, which chops it back down so that -1 < rndm < 1. 

         if(abs(rndm).gt.1.)then
           rndm = amod(rndm, 1.)
         endif

         return
      end

c*******************************************************************

