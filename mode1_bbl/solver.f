C-This file contains all subroutines related to the Navier-Stokes solver.



c*********************************************************************  
      subroutine convec(istep,istart,ispecfilternl)
c*********************************************************************  
c
      include 'dim.h'

!      parameter (dtanim = 5.0) ! How often is data dumped out
c
      include 'comflow.h'
      include 'commatrix.h'
      include 'comtime.h'
      include 'comwave.h'
      include 'comleng.h'
c      include 'comparamphys.h'
      include 'comcoeff.h'
      include 'comparamphys.h'
      include 'comfft.h'
      include 'commass.h'
      include 'comscratch.h'
      include 'comgridptxyz.h'
      include 'commeanv.h'
      include 'comrotpressbc.h'
      include 'comsubd.h'
      include 'comtimeadvance.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'
      include 'dimgrad.h'
      include 'complexflow.h'
      include 'complexgrad.h'
c
      dimension uwork(5,nxpp,ny,nz) ! ,tb(nxpp,ny,nz)

      dimension xdudz(nxpp,ny,nz),xdwdx(nxpp,ny,nz),xdtdz(nxpp,ny,nz)

      complex us(nzloc),vs(nzloc),ws(nzloc)
      complex usn(nzloc),vsn(nzloc),wsn(nzloc)
      complex usnm(nzloc),vsnm(nzloc),wsnm(nzloc)
      complex rot,rotn,rotnm,rotf,rx
      character*80 fout
      character*9 fmt2
      character*30 fspec
     
   
c
      include 'equivsuvwt.h'
      include 'equivscratch.h'
      include 'equivflow.h'
c     include 'equivflowprev.h'
      include 'equivgrad.h'
      common /Dump/ Tdump,Tx
      common /anim/ianim,ianimc      
      common /deriv/xdudz,xdwdx,xdtdz
! JMS -common
      common /dtanimation/ dtanim,fwave 


      rx = -1.0*cmplx(0.,1.)      
c
c set mean w to zero
c set mean v to zero (**PD: Test for wake meandering**)
c
      do k=1,nz
c      vc(1,1,k)=(0.,0.)
       wc(1,1,k)=(0.,0.)
      end do
c


      call copy(divg,ox,ntot)
      call dpdz(u,dudz)
      call dpdz(v,dvdz)
      call dpdz(w,dwdz)
      call dpdz(temp,dtdz)
       
c     do k=1,nz
c       write(*,*) u(24,48,k)
c       write(*,*) z(k),u(1,1,k),dudz(1,1,k) ! ,dvdz(1,1,k),dwdz(1,1,k)
c     enddo

c     write(*,*) '///////////'
c
      do k = 1,nz                                                   
        tempmean(k) = temp(1,1,k)
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

c        if ((i.eq.2).and.(j.eq.2)) write(*,*) dudxc(i,j,k)

        enddo
       enddo
      enddo
c
c  spectral-to-physical
c
      call  horfft (u,1,nz)
      call  horfft (v,1,nz)
      call  horfft (w,1,nz)
      call  horfft (dudx,1,nz)
      call  horfft (dudy,1,nz)
      call  horfft (dudz,1,nz)
      call  horfft (dvdx,1,nz)
      call  horfft (dvdy,1,nz)
      call  horfft (dvdz,1,nz)
      call  horfft (dwdx,1,nz)
      call  horfft (dwdy,1,nz)
      call  horfft (dwdz,1,nz)
      call  horfft (ox,1,nz)
      call  horfft (temp,1,nz)
      call  horfft (dtdx,1,nz)
      call  horfft (dtdy,1,nz)
      call  horfft (dtdz,1,nz)

      call zero(su,ntot)
      call zero(sv,ntot)
      call zero(sw,ntot)
      call zero(st,ntot)

c************************************
                xdudz=dudz
                xdwdx=dwdx          
                xdtdz=dtdz          

C************************************

c-Dump out data for animations
c-Is integer of current time a multiple of dtanim     
c         open(987,file='testvort')
c      if (mod( float(int(t)),dtanim).eq.0) then 
         if ((t-Tdump).ge.dtanim) then
           Tdump=t
c        if (ianim.eq.0) then
c          ianim = 1
c************ Changing file extensions if Restarting the solution
c*********** Correct extension is 1+int(t/dtanim)
C*********** A.M.A. June-04-2008 
c        if (istart.eq.1) then
c          ianimc =1+int(t/dtanim)
c        else
C************* Extension if there is no Restart ************
         ianimc=ianimc+1
c        endif

c         timeout = t
          if(istart.eq.1) write(*,*)'ianimc=',ianimc
          timeout = ianimc
          fmt2='(i5.5,a4)'
          if (nint(timeout).ge.10)  fmt2='(i5.5,a4)'
          if (nint(timeout).ge.100) fmt2='(i5.5,a4)'
          if (nint(timeout).ge.1000) fmt2='(i5.5,a4)'
          if (nint(timeout).ge.10000) fmt2='(i5.5,a4)'
          fspec = 'vortanimout_'
          fout = fspec
          write(unit=fout(13:23),fmt=fmt2) ianimc,'.dat'
          
          write(*,*) 'Dumping out animation file to',
     >    fout,'at',t
          open(500,file=fout)  
          write(500,*) 'ZONE F=POINT, I=',nx,', J=',nz
           espi=1.e-20  !a very small number to avoid inf Rich nu
           dergrad= rhograd(zlen,zlen,grav,rho0,brunt)
          do ks=1,nsubd
            do kloc=1,nzloc
              k = (ks-1)*nzloc + kloc
              zz = z(k)
              shearderv = uo*fcshear*dfsheardz(zz,zdim,zsh,xacu)
              do i=1,nx
                xx = x(i)
    
!                drhomean = abs( rhomean(z(nz),zlen,grav,rho0,brunt) -
!     >               rhomean(z(1),zlen,grav,rho0,brunt) )
 
                vorty = dwdx(i,j,k) - dudz(i,j,k)     

                write(500,'(1x,6(F25.16,3x))') xx,zz,u(i,j,k),w(i,j,k),
     >         temp(i,j,k),
!    >         (-grav*dtdz(i,j,k)/rho0),
     >         temp(i,j,k)+rhomean(zz,zlen,grav,rho0,brunt)
!     >      max((dudz(i,j,k)+shearderv)**2,epsi)

               enddo
             end do
           end do

!            write(*,*)'CL.Rich No.',   
!     > brunt**2/(uo*fcshear*dfsheardz(0.79,zdim,zsh,xacu))**2  ! good!

           close(500)

           open(67,file='animdump.log',access='append')
           write(67,*) 'Animation data dumped out at time',t,
     >              'vortanimout_',ianimc,'.dat','at step=',istep
           close(67)

c-Special case where dtanim = dt
c           if (dtanim.eq.dt) ianim=0

c       else
c-File has been dumped out already within this interval 
c        continue 

c        endif
c      else
c-If integer of given time is not a multiple of dtanim but the
c-animation flag is still on this means that we have progressed
c-to the integer time immediately following one where an animation
c-dump was performed. Set the animation flag to 0
c        if (ianim.eq.1) ianim=0
c      continue

      endif

C-Use non-conservative form suggested by Dan Chan
C-30+ subsequent lines from D. Bogucki's code have been eliminated
c
c
C-Employ explicitly conservative form for momentum equations
C-Adapted from LES code
C-PD: 2/14/02. DTDz never used in convec. Also got rid of eddyvis variable.
      call rhsfm(dudx,dudy,dudz,dvdx,dvdy,dvdz, 
     +           dwdx,dwdy,dwdz,dtdx,dtdy,dtdz,
     +           suc,svc,swc,stc,tc,istep)

c
C-Now calculate volume and surface integrals of various rhs expressions
C-First calculate volume integrals
c     call numintspec(alpha,beta,volintu,volintv,volintw)
c
c  [u] T
c
c

C-Save the advective term estimates for timestep (n)
C-in temporary arrays. Do same for velocity/temperature
C-values at (n)

      do k=1,nz  
        do j=1,ny
          do i=1,nxpp

            uwork(1,i,j,k) = su(i,j,k)
            uwork(2,i,j,k) = sv(i,j,k)
            uwork(3,i,j,k) = sw(i,j,k)
            uwork(4,i,j,k) = st(i,j,k)
            uwork(5,i,j,k) = u(i,j,k)
            ox(i,j,k) = v(i,j,k)
            oy(i,j,k) = w(i,j,k)
            oz(i,j,k) = temp(i,j,k)

c           if ((i.eq.1).and.(j.eq.1)) 
c    >      write (*,*) k,su(i,j,k),sv(i,j,k),sw(i,j,k)
          enddo
        enddo
       enddo

C-Compute everything necessary for improved pressure BC's.
C-Calculate NON-LINEAR term contribution of w-velocity to vertical
C-pressure BC.
C-Bottom boundary
      do j=1,ny
        do i=1,nxpp

          wnl(i,j,1) =  fac1*sw(i,j,1)
     >                 + fac2*swn(i,j,1)
     >                 + fac3*swnm(i,j,1)
C-Top Boundary
          wnl(i,j,2) = fac1*sw(i,j,nz)
     >                 + fac2*swn(i,j,nz)
     >                 + fac3*swnm(i,j,nz)

c         write(*,*) i,j,su(i,j,nz),wnl(i,j,2)
        enddo
      enddo

c
c----------------------------------------------------------------            
c
c-Now express the RHS of the explicit step (NL-term advancement)
c-as the AB3 approximation of the advective term. 
c
c----------------------------------------------------------------             
c
      do k=1,nz                                                    
       do j=1,ny                                                     
        do i=1,nxpp                                                    

         su(i,j,k) =  fac1*su(i,j,k) 
     >              + fac2*sun(i,j,k)
     >              + fac3*sunm(i,j,k)  

         sv(i,j,k) =   fac1*sv(i,j,k) 
     >              + fac2*svn(i,j,k) 
     >              + fac3*svnm(i,j,k) 

         sw(i,j,k) =   fac1*sw(i,j,k) 
     >              + fac2*swn(i,j,k)
     >              + fac3*swnm(i,j,k) 

         st(i,j,k) =   fac1*st(i,j,k) 
     >              + fac2*stn(i,j,k)
     >              + fac3*stnm(i,j,k) 

        enddo
       enddo
      enddo

C-Now focus on bottom and top elements to calculate rotational
C-component of VISCOUS term in pressure BC's
      do j=1,ny
        ii=0
        do i=1,nxpp,2
          ii = ii+1

          xwin =  xw(ii)
          ywin =  yw(j)
C-BOTTOM subdomain
          ks = 1
          ybar = zh(ks)
          do kloc = 1,nzloc
            k = (ks-1)*nzloc + kloc

            us(kloc) = cmplx ( u(i,j,k),u(i+1,j,k) )
            vs(kloc) = cmplx ( v(i,j,k),v(i+1,j,k) )
            ws(kloc) = cmplx ( w(i,j,k),w(i+1,j,k) )

            usn(kloc) = cmplx ( un(i,j,k),un(i+1,j,k) )
            vsn(kloc) = cmplx ( vn(i,j,k),vn(i+1,j,k) )
            wsn(kloc) = cmplx ( wn(i,j,k),wn(i+1,j,k) )

            usnm(kloc) = cmplx ( unm(i,j,k),unm(i+1,j,k) )
            vsnm(kloc) = cmplx ( vnm(i,j,k),vnm(i+1,j,k) )
            wsnm(kloc) = cmplx ( wnm(i,j,k),wnm(i+1,j,k) )

          enddo

          call rotvisc(us,vs,ws,xwin,ywin,ybar,rot,1)
          call rotvisc(usn,vsn,wsn,xwin,ywin,ybar,rotn,1)
          call rotvisc(usnm,vsnm,wsnm,
     >                 xwin,ywin,ybar,rotnm,1)

          rotf = fac1*rot + fac2*rotn + fac3*rotnm

          wviscrot(i,j,1) = -xnu*real(rotf)
          wviscrot(i+1,j,1) = -xnu*real(rx*rotf)

C-TOP subdomain
          ks = nsubd
          ybar = zh(ks)
          do kloc = 1,nzloc
            k = (ks-1)*nzloc + kloc

            us(kloc) = cmplx ( u(i,j,k),u(i+1,j,k) )
            vs(kloc) = cmplx ( v(i,j,k),v(i+1,j,k) )
            ws(kloc) = cmplx ( w(i,j,k),w(i+1,j,k) )

            usn(kloc) = cmplx ( un(i,j,k),un(i+1,j,k) )
            vsn(kloc) = cmplx ( vn(i,j,k),vn(i+1,j,k) )
            wsn(kloc) = cmplx ( wn(i,j,k),wn(i+1,j,k) )

            usnm(kloc) = cmplx ( unm(i,j,k),unm(i+1,j,k) )
            vsnm(kloc) = cmplx ( vnm(i,j,k),vnm(i+1,j,k) )
            wsnm(kloc) = cmplx ( wnm(i,j,k),wnm(i+1,j,k) )

          enddo

          call rotvisc(us,vs,ws,xwin,ywin,ybar,rot,nzloc)
          call rotvisc(usn,vsn,wsn,xwin,ywin,ybar,rotn,nzloc)
          call rotvisc(usnm,vsnm,wsnm,
     >                 xwin,ywin,ybar,rotnm,nzloc)

          rotf = fac1*rot + fac2*rotn + fac3*rotnm

          wviscrot(i,j,2) = -xnu*real(rotf)
          wviscrot(i+1,j,2) = -xnu*real(rx*rotf)

!         if ((i==1).and.(j==1))
!    >     write(*,*) wviscrot(i,j,1),wviscrot(i+1,j,1)
!         if ((i==1).and.(j==1))
!    >     write(*,*) wviscrot(i,j,2),wviscrot(i+1,j,2)
c         if ((ii.eq.2).and.(j.eq.2)) write(*,*)
c    >rotf,rot,rotn,rotnm
        enddo
      enddo

C-Final loop: Update u,v,w,T  at intermediate level with AB3 of 
C-n.l. term and BDF of variables at preceding timesteps

! JMS -Test
!      open(54,file='testu.dat')

      do k=1,nz
        do j=1,ny
          do i=1,nxpp            
              xdummy = u(i,j,k)
              u(i,j,k) = dt*su(i,j,k)
     >              + bdf1*u(i,j,k)
     >              + bdf2*un(i,j,k)
     >              + bdf3*unm(i,j,k)
!              if (k.eq.nz) then
!              write(54,'(1x,3(e14.7,3x))') z(k),xdummy, u(i,j,k)
!              endif 
              v(i,j,k) = dt*sv(i,j,k)
     >              + bdf1*v(i,j,k)
     >              + bdf2*vn(i,j,k)
     >              + bdf3*vnm(i,j,k)

              w(i,j,k) =  dt*sw(i,j,k)
     >              + bdf1*w(i,j,k)
     >              + bdf2*wn(i,j,k)
     >              + bdf3*wnm(i,j,k)

              temp(i,j,k) = dt*st(i,j,k)
     >              + bdf1*temp(i,j,k)
     >              + bdf2*tn(i,j,k)
     >              + bdf3*tnm(i,j,k)
            
          enddo
        enddo
      enddo 
c-------------------------------------------------------------------
c-NOTE: As a warning, the "subscripts" used for the various advective
c-terms may be a bit confusing. The example below, may help:
c-u(...) = velocity at level (n)
c-un(...) = velocity at level (n-1)
c-unm(...) = velocity at level (n-2)
c---------------------------------------------------------------------


c
c  update the convective terms at n-2 and n-1 levels
c  *work is what was originally s* and now is transferred to
c  level n-1
c  *n is transferred to *nm (i.e. n-1 contribution becomes
c  n-2 for next step)
c
        do k=1,nz
          do j=1,ny
            do i=1,nxpp

             
C-Relegate non-linear term at level (n-1) to (n-2)
              sunm(i,j,k) = sun(i,j,k)
              svnm(i,j,k) = svn(i,j,k)
              swnm(i,j,k) = swn(i,j,k)
              stnm(i,j,k) = stn(i,j,k)

C-Relegate non-linear term at level (n) to (n-1)
              sun(i,j,k) = uwork(1,i,j,k)
              svn(i,j,k) = uwork(2,i,j,k)
              swn(i,j,k) = uwork(3,i,j,k)
              stn(i,j,k) = uwork(4,i,j,k)

C-Relegate velocity/temperature at level (n-1) to (n-2)
              unm(i,j,k) = un(i,j,k)
              vnm(i,j,k) = vn(i,j,k)
              wnm(i,j,k) = wn(i,j,k)
              tnm(i,j,k) = tn(i,j,k)

C-Relegate velocity/temperature at level (n) to (n-1)
              un(i,j,k) = uwork(5,i,j,k)
              vn(i,j,k) = ox(i,j,k)
              wn(i,j,k) = oy(i,j,k)
              tn(i,j,k) = oz(i,j,k)

            enddo
          enddo
        enddo             

c-Now if desired call routine that filters data in
c-Fourier and Legendre space after N-L term advancement
        ifourier = 1
        ilegen = 1
        if  (ispecfilternl.eq.1) then
          write(*,*) 'Filtering NON-LINEAR terms'
          call filter_3dlegenfourier(u,ifourier,ilegen,0)
          call filter_3dlegenfourier(v,ifourier,ilegen,0)
          call filter_3dlegenfourier(w,ifourier,ilegen,0)
          call filter_3dlegenfourier(temp,ifourier,ilegen,0)
        endif

c       do k=1,nz
c         do j=1,ny
c           ii = 0 
c             do i=1,nxpp,2
c               ii = ii +1

c                 if ((ii.eq.1).and.(j.eq.1)) then
c                   w(i,j,k) = 0.
c                 else
c                   w(i,j,k) = w(i,j,k) + dt*ra*pr*temp(i,j,k)
c                   w(i+1,j,k) = w(i+1,j,k) + dt*ra*pr*temp(i+1,j,k)
c                   if ((ii.eq.2).and.(j.eq.2)) write(*,*) 
c    >              z(k),ra*pr*temp(i,j,k),ra*pr*temp(i+1,j,k)
c                 endif
c             enddo
c         enddo
c       enddo
           
c


      return
      end             







c****************************************************************
      subroutine rhsfm(dudx,dudy,dudz,dvdx,dvdy,dvdz, ! 
     +                dwdx,dwdy,dwdz,dtdx,dtdy,dtdz,
     +                suc,svc,swc,stc,tc,istep)
c****************************************************************
C-This is a subroutine which implements calculation of nonlinear
C-term in skew-symmetric (conservative form). A variant of what
C-existed in Shari Kimmel's LES code. Lots of Shari's stuff
C-has been removed (**PD: October 2001**)
      include 'dim.h'

      parameter(tgrad=1.)
c
      include 'comflow.h'
      include 'commatrix.h'
      include 'comtime.h'
      include 'comwave.h'
      include 'comparamphys.h'
      include 'comleng.h'
      include 'comscratch.h'
      include 'commass.h'
      include 'dimgrad.h'
c     include 'dimtempgrad.h'
      include 'comsubd.h'
      include 'commeanv.h'
      include 'comgridptxyz.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'
c     character*9 testfile
c     character*5 f1
c-------------------------------------------------------------------


      real skewsym1(nxpp,ny,nz),skewsym2(nxpp,ny,nz)         

C-PD: 2/14/02. These common blocks are commented out because
C-they are only used in complete version of LES estimation model.
c     common /flowf/   ut(ntot),vt(ntot),wt(ntot)
c     common /param/   pr,ra,prt,Cs,delta(nz),ybar,re
c     common /tempst/  tempsw,cmt,ctt  !CHANNEL: Temperature switch and c
c     common /control/ iflag,iprt,itype,imodel,mprt
c          common  /forcef/  phiz(nz)


      complex o1(nxhp,ny,nz),o2(nxhp,ny,nz),o3(nxhp,ny,nz),
     1        tc(nxhp,ny,nz),stc(nxhp,ny,nz),
     4        suc(nxhp,ny,nz),svc(nxhp,ny,nz),swc(nxhp,ny,nz)
      complex skewsym1c(nxhp,ny,nz),skewsym2c(nxhp,ny,nz)     
C-PD: 2/14/02. The dTdx, dTdy don't seem to be used anywhere, so I
C-commented them out.
c     complex dtdxc(nxhp,ny,nz),dtdyc(nxhp,ny,nz)



      equivalence (skewsym1,skewsym1c),(skewsym2,skewsym2c)
      include 'equivscratch.h'
c     equivalence (dtdx,dtdxc),(dtdy,dtdyc)

! JMS -common
      common /dtanimation/ dtanim,fwave







C-**PD: This is a fossil comment from Shari's code** 
c ------ get (u.grad)u by using conservation form
c calculate additional terms for LES--first terms that need d/dz
c include terms requiring d/dz from the second part of the nonlinear term
c the nonlinear term is expressed as
c (u.grad)u= 0.5(u.grad)u + 0.5*grad(u.u) in physical space
c example for u-component:
c     0.5*{d(uu)/dx+d(uv)/dy+d(uw)/dz+
c          u*[d(u+t11)/dx]+v*[d(u+t12)/dy]+w*[d(u+t12)/dz]}




C-**PD: itag is evidently a flag which indicates whether conservative
C-form is employed** I left this for no specific reason :)
      itag = 1 
      if(itag .eq. 1) then


           if(istep.eq.1) then
                open(660,file='phizcutx.dat')
                write(660,*) 'ZONE F=POINT, I=',nx,', J=',nz
                open(111,file='shearprof.dat') 
                write(111,*)'VARIABLES = "shear","z",
     >                            "dsheardz","d2sheardz2"'
           endif
 



c **PD**: The process consists of the three steps:
c a) Set up array tempz in physical space
c b) Employing differentation properties of Legendre
c    polynomials calculate ddz for 
c Update corresponding su, sv or sw term all in one shot
c
c duwdz term in x-momentum equation
      
      do k=1,nz
          do j=1,ny
            do i=1,nxpp

            skewsym1(i,j,k) = -0.5*u(i,j,k)*w(i,j,k)

c           if ((i.eq.1).and.(j.eq.1)) write(*,*) 
c    >      u(1,1,k),w(1,1,k),skewsym1(1,1,k)
          enddo
        enddo
       enddo 

      call dpdz(skewsym1,su)       

c     write(*,*) '////////////////' 
c dvwdz term in y-momentum equation
      do k=1,nz
          do j=1,ny
            do i=1,nxpp

            skewsym1(i,j,k) = -0.5*v(i,j,k)*w(i,j,k)

c           if ((i.eq.1).and.(j.eq.1)) write(*,*) su(1,1,k)

          enddo
        enddo
       enddo

      call dpdz(skewsym1,sv)
    
c dw^2dz term in z-momentum equation
       do k=1,nz
          do j=1,ny
            do i=1,nxpp

            skewsym1(i,j,k) = -0.5*w(i,j,k)*w(i,j,k)

          enddo
        enddo
       enddo

      call dpdz(skewsym1,sw)

c-now update su,sv & sw terms all in one shot    
      do k=1,nz
        zz = z(k)
        do j=1,ny
          do i=1,nxpp

            su(i,j,k)=su(i,j,k)
     >                -0.5*(u(i,j,k)*dudx(i,j,k)
     >                +v(i,j,k)*(dudy(i,j,k)-4.*omz)
     >                +w(i,j,k)*dudz(i,j,k))

            sv(i,j,k)=sv(i,j,k)
     >                -0.5*(u(i,j,k)*(dvdx(i,j,k)+4.*omz)
     >                +v(i,j,k)*dvdy(i,j,k)
     >                +w(i,j,k)*dvdz(i,j,k))

            sw(i,j,k)=sw(i,j,k)
     >                -0.5*(u(i,j,k)*dwdx(i,j,k)
     >                +v(i,j,k)*dwdy(i,j,k)
     >                +w(i,j,k)*dwdz(i,j,k)) 

C-Add GRAVITY term in physical space
c           tfluct = (temp(i,j,k) - tempmean(k))
            sw(i,j,k) = sw(i,j,k) + 
     >                 ((bruntz(zz,zlen,brunt)**2.)/
     >                  rhograd(zz,zlen,grav,rho0,brunt))*
     >                  temp(i,j,k)
           
c           if ((i.eq.nxhp).and.(j.eq.nyh))  write(*,*) z(k),tfluct
c    >write(*,*) su(i,j,k),sv(i,j,k),sw(i,j,k)
          enddo
        enddo
       enddo

c-Add Non-Linear PENALTY terms
c-Contribution appears only at boundaries and interfaces
c-Scan over subdomains in vertical

      omega = 2./((nzloc+1)*(nzloc+2))

      do i=1,nxpp
        do j=1,ny
          do ks=1,nsubd


              ybar = zh(ks)
              kb = (ks-1)*nzloc + 1
              kt = (ks-1)*nzloc + nzloc

C-Determine if I have inflow/outflow at interface.
C-At the walls (w=0) consider this to be an "outflow".
C-Note that I consider the total vertical velocity (perturbation
C- + mean wave) for the soliton problem.
              wb = w(i,j,kb) + wmwv(i,kb)
              wt = w(i,j,kt) + wmwv(i,kt)
              if (kb.eq.1) wb = 0.
              if (kt.eq.nz) wt = 0.
C-Set penalty parameters
C-Lower interface
              if (wb.gt.0.) then
                alpha0 = wb
                tau1 = 0.5/omega*(2./ybar)
              else
                alpha0 = 0.
                tau1 = 0.
              endif

c             if (alpha0.eq.0.) then
c               write(*,*) 'tau1',tau1,alpha0,kb,kt
c             endif

C-Upper interface
              if (wt.lt.0.) then
                gamma0 = abs(wt)
                tau2 = 0.5/omega*(2./ybar)
              else
                gamma0 = 0.
                tau2 = 0.
              endif

c             if (gamma0.eq.0.) then
c               write(*,*) 'tau2',tau2,gamma0,kb,kt
c             endif
c             write(*,*) i,j,kb,kt
c             write(*,*) wb,wt,alpha,gamma,tau1,tau2

C-Now add penalty terms to bottom and top interfaces:
C-U-VEL

              if (kb.eq.1) then
                pcb = 0.
              else
                pcb = alpha0*u(i,j,kb-1)
              endif

              su(i,j,kb) = su(i,j,kb) 
     >      - tau1*(alpha0*u(i,j,kb) - pcb)


              if (kt.eq.nz) then
                pct = 0.
              else
                pct = gamma0*u(i,j,kt+1)
              endif

              su(i,j,kt) = su(i,j,kt) 
     >      - tau2*(gamma0*u(i,j,kt) - pct)

C-V-VEL
              if (kb.eq.1) then
                pcb = 0.
              else
                pcb = alpha0*v(i,j,kb-1)
              endif

              sv(i,j,kb) = sv(i,j,kb)
     >      - tau1*(alpha0*v(i,j,kb) - pcb)


              if (kt.eq.nz) then
                pct = 0.
              else
                pct = gamma0*v(i,j,kt+1)
              endif

              sv(i,j,kt) = sv(i,j,kt)
     >      - tau2*(gamma0*v(i,j,kt) - pct)
C-W-VEL
              if (kb.eq.1) then
                pcb = 0.
              else
                pcb = alpha0*w(i,j,kb-1)
              endif

              sw(i,j,kb) = sw(i,j,kb)
     >      - tau1*(alpha0*w(i,j,kb) - pcb)

              if (kt.eq.nz) then
                pct = 0.
              else
                pct = gamma0*w(i,j,kt+1)
              endif

              sw(i,j,kt) = sw(i,j,kt)
     >      - tau2*(gamma0*w(i,j,kt) - pct)



C-TEMPERATURE
              if (kb.eq.1) then
                pcb = 0.
              else
                pcb = alpha0*temp(i,j,kb-1)
              endif

              st(i,j,kb) = st(i,j,kb)
     >      - tau1*(alpha0*temp(i,j,kb) - pcb)

              if (kt.eq.nz) then
                pct = 0.
              else
                pct = gamma0*temp(i,j,kt+1)
              endif

              st(i,j,kt) = st(i,j,kt)
     >      - tau2*(gamma0*temp(i,j,kt) - pct)
 
          enddo
        enddo
       enddo

       
C-FORCING TERMS DUE TO SOLITON
C-PD: 2/24/04 = By mean here, we mean wave contributions. 

! JMS-test
!         open(53,file='testsu.dat')
        
       do k=1,nz ! Bottom point overlooked because of solid wall
                 ! where xorcing has no meaning

         zz = z(k)
c         xphiz=phiz(zz)

         do j=1,ny
          do i=1,nxpp
             xx = x(i)                     
C-1) Advection of perturbation by mean flow (and free stream/shear)
c-(u,v,w,rho) eqns.
c-CAREFUL: Sign of umax (upstream/downstream propagation)
C- note that dumdx,dumdz,dwmdx,dwmdz,drhomdz,drhomdx have all
C-been set to zero in wave setup.f // A.M.A.
C- Here we are in lab coordinates (no moving frame)A.M.A. 03-11-08
c       write(*,*)'dudz=',dudz(i,j,k),'dwdx',dwdx(i,j,k)            
c           if (t.lt.600)fcshear=0
c           if(t.ge.0)fcshear=1
c            if(i.eq.32)write(*,*)'dtdx=',dtdx(i,j,k),'dwdx',dwdx(i,j,k)
c             fcshear=1.0 
C          write(*,*)'uo',uo,'fcshear',xomega

          advmean = umwv(i,k) +
     >           uo*fcshear*fshear(zz,zdim,zsh,xacu)
!            if (zz.eq.z(nz)) then  
!             write(53,'(1x,2(e14.7,3x))') zz,su(i,j,k)
!            endif 

! trying to find why u is getting set to zero at top boundary 

             xdummy = su(i,j,k)
             su(i,j,k) = su(i,j,k) - advmean*dudx(i,j,k) 
     >                             - wmwv(i,k)*dudz(i,j,k)
!             if (zz.eq.z(nz)) then 
!             write(53,'(1x,3(e14.7,3x))') zz,xdummy,su(i,j,k)
!             endif 
c**********this accounts for -U(z)du'/dx source term
             sv(i,j,k) = sv(i,j,k) - advmean*dvdx(i,j,k) 
     >                             - wmwv(i,k)*dvdz(i,j,k)

             sw(i,j,k) = sw(i,j,k) - advmean*dwdx(i,j,k) 
     >                             - wmwv(i,k)*dwdz(i,j,k)
c************** this accounts for -U(z)dw'/dx term

             st(i,j,k) = st(i,j,k) - advmean*dtdx(i,j,k)
     >                             - wmwv(i,k)*dtdz(i,j,k)
c*************** this accounts for -U(z)drho'dx term

         

C-2) Advection of mean flow (take shear into account) by perturbation
c-(u,w,rho) eqns.
             
             dusheardz = uo*fcshear*dfsheardz(zz,zdim,zsh,xacu)
           
             su(i,j,k) = su(i,j,k) - u(i,j,k)*dumdx(i,k)
     >                             - w(i,j,k)*dumdz(i,k)
     >                             - w(i,j,k)*dusheardz 

             if (zz.eq.z(nz)) then 
             write(53,'(1x,3(e14.7,3x))') zz,xdummy,su(i,j,k)
             endif 

c******** this accounts for -w'dU(z)/dz term
             sw(i,j,k) = sw(i,j,k) - u(i,j,k)*dwmdx(i,k)
     >                             - w(i,j,k)*dwmdz(i,k)

             st(i,j,k) = st(i,j,k) - u(i,j,k)*drhomdx(i,k)
     >                             - w(i,j,k)*drhomdz(i,k)


c        write(*,*)'u0',uo,'wmega',xomega,xkcrit
c- visouc contribution of the shear flow A.M.A. 03-05-2008
c%%% The viscous term leads to inconsistent system
c if the wave forcing was turned off. IF the forcing is turned
c on this term leads to horizontal disturbances around the center
c of the current. A.M.A. 06-17-2008 
c         d2usheardz = uo*fcshear*d2fsheardz(zz,zdim,zsh,xacu)
            
c         su(i,j,k) = su(i,j,k) + xnu*d2usheardz
c******** this accounts for nu d2Udz2 term
               
c remember d2rhomdz2=0.0  //A.M.A.
        st(i,j,k) = st(i,j,k) + 
     >        xkappa*d2rhomeandz2(zz,zlen,grav,rho0,brunt)

       
C-VIRTUAL PADDLE !!!
C-8) Forcing due to periodically oscillating localized source.
C-   This source emits a wave-packet downwards (See Slinn & Riley, JCP98)
c           open(679,file='testsolver1.dat')
           

           forceu_rhs = forceu(xx,zz,t,
     >                  xacrit,xkcrit,xmcrit,xomega,
     >                  zdim,umax,fwave)
           forcew_rhs = forcew(xx,zz,t,
     >                  xacrit,xkcrit,xmcrit,xomega,
     >                  zdim,umax,fwave)
           forcerho_rhs = forcerho(xx,zz,t,
     >                  xacrit,xkcrit,xmcrit,xomega,
     >                  zdim,umax,
     >                  brunt,rho0,grav,fwave)

           su(i,j,k) = su(i,j,k) + forceu_rhs

           sw(i,j,k) = sw(i,j,k) + forcew_rhs

           st(i,j,k) = st(i,j,k) + forcerho_rhs

c           write(768,*)'umax,forceu_rhs,
c     >   forcew_rhs,forcerho_rhs,xkappa',umax,forceu_rhs,forcew_rhs,
c     >   forcerho_rhs,xkappa  
c            open(565,file='testenvelope')
c           if((i==10).and.(k==20)) then
c           write(565,*)t,phit(t)
c           endif

           xphiz = phiz(xx,zz,xkcrit,xmcrit,zcen,zdim,xbcrit,xlen)
     >                  *bruntz(zz,zlen,brunt)/0.39
           xphizmin =   phiz(xx,zz,xkcrit,xmcrit,zcen,zdim,xbcrit,xlen)
       
           if(xphiz.le.xphizmin) then
              xphiz = xphizmin
           endif 

           if(xphiz.gt.4.) then
              xphiz = 4.
           endif 


           if((istep.eq.1).and.(i.le.nx)) then 
                write(660,'(3(E12.6,2x),I)') xx,zz,xphiz,i
                if(i.eq.1) then
                write(111,*) uo*fshear(zz,zdim,zsh,xacu),zz,
     > uo*dfsheardz(zz,zdim,zsh,xacu),uo*d2fsheardz(zz,zdim,zsh,xacu)
                endif
           endif

!!
!-PD-Lyon, 7/1/15: De-activated all sponge layer functions.
CCC**************Rayleigh absorber****************************
!          su(i,j,k) = su(i,j,k) - u(i,j,k)*xphiz
 
!          sw(i,j,k) = sw(i,j,k) - w(i,j,k)*xphiz

!          st(i,j,k) = st(i,j,k) - temp(i,j,k)*xphiz
c*************************************************************

c          open(564,file='testabsorber.dat')
c          if(i.eq.32)write(564,*)zz,xphiz

           enddo
         enddo
         
       enddo

       write(*,*) 'Value of FORCING function ',phit(t,xomega,fwave)

       if(istep.eq.1) then
            close(660)
            close(111)
       endif


c *** terms in x and y direction
c The steps are
c a) In physical space calculate the contribution
c    to each velocity component equation
c b) Perform fft to convert to spectral space
c c) Differentiate and update corresponding suc,svc and swc term
c
c

c u-momentum equation
c du^2dx and duvdy terms
c define
      do k=1,nz
        do j=1,ny
          do i=1,nxpp
            skewsym1(i,j,k) = -0.5*u(i,j,k)*u(i,j,k)
            skewsym2(i,j,k) = -0.5*u(i,j,k)*v(i,j,k)
          enddo
        enddo
      enddo

c physical to spectral space to calculate ddx and ddy
      call  horfft (su,-1,nz)
      call  horfft (skewsym1,-1,nz)
      call  horfft (skewsym2,-1,nz)

c Now differentiate and update su (suc in spectral space)
      do k=1,nz
       do j=1,ny
        do i=1,nxhp      
c         if ((i.eq.1).and.(j.eq.1)) write(*,*) suc(i,j,k)
         suc(i,j,k) = suc(i,j,k)+(0.,1.)*
     +        (xw(i)*skewsym1c(i,j,k)+yw(j)*skewsym2c(i,j,k))
c        write(100) suc(i,j,k)
        enddo
       enddo
      enddo

c v-momentum equation
c duvdy and dv^2dy terms
c define
      do k=1,nz
        do j=1,ny
          do i=1,nxpp
            skewsym1(i,j,k) = -0.5*u(i,j,k)*v(i,j,k)
            skewsym2(i,j,k) = -0.5*v(i,j,k)*v(i,j,k)
           enddo
        enddo
      enddo

c physical to spectral space to calculate ddx and ddy
      call  horfft (sv,-1,nz)
      call  horfft (skewsym1,-1,nz)
      call  horfft (skewsym2,-1,nz)

c Now differentiate and update sv (svc in spectral space)
      do k=1,nz
       do j=1,ny
        do i=1,nxhp
         svc(i,j,k) = svc(i,j,k)+(0.,1.)*
     +        (xw(i)*skewsym1c(i,j,k)+yw(j)*skewsym2c(i,j,k))
c        write(100) svc(i,j,k)
        enddo
       enddo
      enddo

c z-momentum equation
c duwdy and dvwdy terms
c define
      do k=1,nz
        do j=1,ny
          do i=1,nxpp
            skewsym1(i,j,k) = -0.5*u(i,j,k)*w(i,j,k)
            skewsym2(i,j,k) = -0.5*v(i,j,k)*w(i,j,k)
          enddo
        enddo
      enddo

c physical to spectral space to calculate ddx and ddy
      call  horfft (sw,-1,nz)
      call  horfft (skewsym1,-1,nz)
      call  horfft (skewsym2,-1,nz)

c Now differentiate and update sw (swc in spectral space)
      do k=1,nz
       do j=1,ny
        do i=1,nxhp
         swc(i,j,k) = swc(i,j,k)+(0.,1.)*
     +        (xw(i)*skewsym1c(i,j,k)+yw(j)*skewsym2c(i,j,k))
c        write(100) swc(i,j,k)
        enddo
       enddo
      enddo
  
      else

        write(*,*) 'Do not go gentle into that good night !'

      end if


c ------ temperature terms

C-PD: 5/28/03: We are now solving for the temperature PERTURBATION
c
c  [u] T + q
c

c-Calculate t*
c
c-d(wT)/dz term in temperature equation
      do k=1,nz
        do j=1,ny
          do i=1,nxpp
            skewsym1(i,j,k)=-w(i,j,k)*temp(i,j,k) 
c    >                      -w(i,j,k)*z(k)/zlen
           enddo
        enddo
      enddo

C-Careful (PD: 5/27/03). There was a major error here.
C-st (already defined was overwritten) and then set equal
C-to d(wT)/dz + wT (???)
      call dpdz(skewsym1,ox)

c-update st term
      do k=1,nz
        do j=1,ny
          do i=1,nxpp

            zz=z(k) 
c           st(i,j,k)=st(i,j,k)+skewsym1(i,j,k)
            st(i,j,k)=st(i,j,k)+ox(i,j,k)  - 
     >                w(i,j,k)*rhograd(zz,zlen,grav,rho0,brunt)
          enddo
        enddo
      enddo

c-Terms d(uT)/dx and d(vT)/dy in temperature equation
c
      do k=1,nz
        do j=1,ny
          do i=1,nxpp
            skewsym1(i,j,k)=-u(i,j,k)*temp(i,j,k)
            skewsym2(i,j,k)=-v(i,j,k)*temp(i,j,k)
          enddo
        enddo
      enddo

c-convert st and the two above terms to spectral space 
      call horfft (st,-1,nz)
      call horfft (skewsym1,-1,nz)
      call horfft (skewsym2,-1,nz)

c-convert temperature to Fourier space
      call  horfft (temp,-1,nz)

c-Now differentiate and update st (stc in spectral space)
c-Add Boussinesq term to rhs of w-velocity equation
      do k=1,nz
       do j=1,ny
        do i=1,nxhp
         stc(i,j,k) = stc(i,j,k)+(0.,1.)*
     +        (xw(i)*skewsym1c(i,j,k)+yw(j)*skewsym2c(i,j,k))
C-Comment out the following line, as we now work with temperature perturbations
C-on RHS
C-GRAVITY TERM
C-PD:2/24/04: Now set in physical space to accomodate arbitrarily
C-BV freq. profile
c        if ((i.eq.1).and.(j.eq.1)) then
c          swc(i,j,k) = cmplx(0.,0.)
c        endif
c        else
c          swc(i,j,k) = swc(i,j,k) + ra*pr*tc(i,j,k)
c          if ((i.eq.1).and.(j.eq.1)) write(*,*) suc(i,j,k)
c        endif
 
c-PD: 5/26/03. Add very specific statement about (0,0) mode for
C-w velocity
c        if ((i.eq.1).and.(j.eq.1))  then
c          swc(i,j,k) = cmplx(0.,0.)
c        endif

        enddo
       enddo
      enddo      
 

      call horfft (u,-1,nz)
      call horfft (v,-1,nz)
      call horfft (w,-1,nz)
           
      return
      end


C*********************************************************
       subroutine rotvisc(us,vs,ws,xw,yw,ybar,rot,k0)
C*********************************************************
c
c      calculates the vertical component of rot(rot(w)) at
c      bottom & top boundaries
c      k0 = indicates whether final results is desired at
c      at bottom or top point of subdomain
c 
       include 'dim.h'
c
       include 'comleng.h'
       include 'comsubd.h'
       complex us(nzloc),vs(nzloc),ws(nzloc)      
       dimension fui(nzloc),fur(nzloc),dfui(nzloc),dfur(nzloc)
       dimension fvi(nzloc),fvr(nzloc),dfvi(nzloc),dfvr(nzloc)
       complex rot,rx

       rx = -1.0*cmplx(0.,1.)

C-Assign values to temporary arrays to calculate dudz and dvdz
       do kloc = 1,nzloc
         fur(kloc) = real(us(kloc))
         fui(kloc) = real(rx*us(kloc))
         fvr(kloc) = real(vs(kloc))
         fvi(kloc) = real(rx*vs(kloc))
       enddo

C-Calculate vertical derivative in specific subdomain
C-PD: CHANGE THIS SOON: Calculate derivative only at k0 !
       call dpdzcol(fur,dfur,d,zpts,nzm,ybar,1,nzloc)
       call dpdzcol(fui,dfui,d,zpts,nzm,ybar,1,nzloc)

       call dpdzcol(fvr,dfvr,d,zpts,nzm,ybar,1,nzloc)
       call dpdzcol(fvi,dfvi,d,zpts,nzm,ybar,1,nzloc)

C-Compute rotational term

       rot = +cmplx(xw*xw,0)*ws(k0)  
     >       +cmplx(yw*yw,0)*ws(k0) 
     >       +cmplx(0,xw)*cmplx(dfur(k0),dfui(k0))
     >       +cmplx(0,yw)*cmplx(dfvr(k0),dfvi(k0))

       return 
       end

C*********************************************************
       function tref(zz,zlen)
C*********************************************************
C-PD: 5/28/03. Specifies reference profile for temperature
C-to be added to calculated temperature perturbations whenever
C-necessary
       tref = zz/zlen

       return
       end

