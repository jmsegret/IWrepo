c*********************************************************************  
       subroutine zero(a,n)
c*********************************************************************  
       dimension a(*)
       do i=1,n
        a(i) = 0.
       enddo
       return
       end




c     
c*********************************************************************  
      subroutine copy(a,b,n)                                              
c*********************************************************************  
      dimension a(*), b(*)                    
c                                                                       
      do i = 1,n                                                   
       b(i) = a(i)                                               
      enddo
      return                                                            
      end                                                               






c
c-----------------------------------------------------
      function random()
c-----------------------------------------------------
c  Routine returns a pseudo-random number between 0-1. 
c-----------------------------------------------------
      integer m, i, md, seed
c      double precision fmd

c      data m/25173/,i/13849/,md/65536/,fmd/65536.d0/,seed/17/
      data m/25173/,i/13849/,md/65536/,fmd/655360/,seed/17/

      save seed

      seed   = mod(m*seed+i,md)
      random = (10.*seed/fmd -0.5)*2. ! -0.5
c     random = random*0.1
      return
      end




c
c---------------------------------------------------------------
      subroutine nusseltsub
c---------------------------------------------------------------
c
      include 'dim.h'
c
      include 'comflow.h'
      include 'comtime.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'comcoeff.h'
      include 'comfft.h'
      include 'comnusselt.h'


      complex uc(nxhp,ny,nz),vc(nxhp,ny,nz),wc(nxhp,ny,nz),
     1        tc(nxhp,ny,nz)
c     
C**PD: ALERT. May have to change the dimensionalization of this.
      dimension dtdz(nxpp,ny,nz)
      real nuss
c
      equivalence (u,uc),   (v,vc),   (w,wc),   (temp,tc)
c-------------------------------------------------------
c
c flow data u,v,w,temp are in the spectral representation
c
c mean temperature (the zero mode in spectral space)
      
      do i=1,nxpp
        do j=1,ny
          do k=1,nz
            tmean(k)=temp(1,1,k)    
          enddo
        enddo
       enddo

c 
c-----  calculate Nusselt number ----------------------
c
c  go to physical space
      call horfft(temp,1,nz)
      call horfft(w,1,nz)
c
c      call dpdz(temp,dtdz,d,nzm,1.0)
c
      nxyz = nx*ny*nz
c
       kkl=1
       kku=nz
       tavl=0.
       tavu=0.
      do j=1,ny
        do i=1,nx
        tavl=tavl+temp(i,j,kkl)
        tavu=tavu+temp(i,j,kku)
         end do
       end do
       tavl=tavl/float(nx*ny) !average temperature at the lower wall
       tavu=tavu/float(nx*ny) !average temperature at the upper wall 
       delt=abs(tavl-tavu)
c
      do k=1,nz
        nuss(k)=0.
        tder(k)=0.
         do j=1,ny   
            do i=1,nx
              nuss(k)=nuss(k) - dtdz(i,j,k) + w(i,j,k)*temp(i,j,k)
        tder(k)=tder(k) + dtdz(i,j,k)
       enddo
      enddo 
       nuss(k)=nuss(k)/float(nx*ny)/delt
       tder(k)=tder(k)/float(nx*ny) 
      enddo
c
c
c  back to spectral space
      call horfft(temp,-1,nz)
      call horfft(w,-1,nz)
c
      return
      end




      
c---------------------------------------------------------------
      subroutine stressb
c---------------------------------------------------------------
c
      include 'dim.h'
c
      include 'comflow.h'
      include 'comtime.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'comcoeff.h'
      include 'comfft.h'
      include 'comscratch.h'
      include 'comsubd.h'
      include 'comgridptxyz.h'

      complex o1(nxhp,ny,nz),o2(nxhp,ny,nz),o3(nxhp,ny,nz),unit
c     real o1(nxhp,ny,nz),o2(nxhp,ny,nz),o3(nxhp,ny,nz)
      real ur(nxpp,ny,nz),vr(nxpp,ny,nz),wr(nxpp,ny,nz)
      dimension deru(nz),derv(nz),derw(nz)
      dimension ozz(nxpp,ny,nz)
      dimension tmp1r(nzloc),tmp1i(nzloc),tmp2r(nzloc),tmp2i(nzloc)
      dimension dtmp1r(nzloc),dtmp1i(nzloc),dtmp2r(nzloc),dtmp2i(nzloc) 

      equivalence (ox,o1), (oy,o2),(oz,o3)
      equivalence (u,ur),(v,vr),(w,wr)

c-------------------------------------------------------
c
c flow data u,v,w,temp are in the spectral representation
c 
c--- calculate derivatives with respect to z
c
c  go to physical space
      call horfft(u,1,nz)
      call horfft(v,1,nz)
      call horfft(w,1,nz)      
c
      call dpdz(u,ox)
      call dpdz(v,oy)
      call dpdz(w,oz)
   
      call dpdz(oz,ozz)

C-AN alternative method of calculating d2wdz2 at top boundary
      ks=1
      ybar = zh(ks)

      do  kloc=1,nzloc
        k = (ks-1)*(nzloc-1) + kloc
        tmp1r(kloc) = w(nxh,nyh,k)
      enddo

      call dpdzcol(tmp1r,dtmp1r,d,zpts,nzm,ybar,1,nzloc)
      call dpdzcol(dtmp1r,dtmp2r,d,zpts,nzm,ybar,1,nzloc)
c     write(*,*) 'Alternative estimate for d2wdz2 at top',
      write(*,*) 'Alternative estimate for dwdz at bottom',
     >dtmp1r(1)
c 
      print *,'u, v, w at the upper surface'
      write(6,11) ur(nxh,nyh,nz),vr(nxh,nyh,nz),wr(nxh,nyh,nz)
      print *,'u, v, w at the lower surface'
      write(6,11) ur(nxh,nyh,1),vr(nxh,nyh,1),wr(nxh,nyh,1)
      print *,'dudz, dvdz, dwdz at the upper surface'
      write(6,10) ox(nxh,nyh,nz),oy(nxh,nyh,nz),oz(nxh,nyh,nz)
c     write(*,*) ox(nxh,nyh,nz),oy(nxh,nyh,nz)
      print *,'dudz, dvdz, dwdz at the lower surface'
      write(6,10) ox(nxh,nyh,1),oy(nxh,nyh,1),oz(nxh,nyh,1)
      print *, 'd2wdz2 at upper surface',ozz(nxh,nyh,nz)
 10   format(1x,6e12.3)       
 11   format(1x,3e12.3)

c     do j=1,ny
c       write(*,*) ox(nxh,j,nz),oy(nxh,j,nz)
c     enddo

c     do k=1,nz
c       write(*,*) w(nxh,nyh,k),oz(nxh,nyh,k),ozz(nxh,nyh,k)
c     enddo
c-Generate test file with output of contours of dudz & dvdz at top
c-surface
      pi2 = 8.0*atan(1.0)
      xl = pi2/alpha
      yl = pi2/beta

      open(600,file='stressbc.dat',status='unknown')
      write(600,*) 'VARIABLES = "x","y","dudz","dvdz"'
      write(600,*) 'ZONE F=POINT, I=',nx,',J=',ny

      do j=1,ny 
        do i=1,nx
          write(600,'(2x,4(E10.4,1x))')
     >x(i),y(j),ox(i,j,nz),oy(i,j,nz)
c         write(600,*) x(i),y(j),ox(i,j,nz),oy(i,j,nz)
        enddo
      enddo
      
      close(600)
        
c

      go to 35

c     gradients of average velocity        
      
      do k=1,nz
         deru(k)=0.
         derv(k)=0.
         derw(k)=0.
         do j=1,ny   
            do i=1,nx
               derv(k)=derv(k) + ox(i,j,k)
               derv(k)=derv(k) + oy(i,j,k)
               derw(k)=derw(k) + oz(i,j,k)       
            enddo
         enddo 
         deru(k)=deru(k)/float(nx*ny) 
         derv(k)=derv(k)/float(nx*ny) 
         derw(k)=derw(k)/float(nx*ny)        
      enddo
      print *,'dudz, dvdz, dwdz ' 
      
      do k=1,nz
         kr=nz-k+1              
C     reverse index for printing from top to bottom
         zz=0.5*(zpts(kr-1)+1.0)
         write(6,92) kr,zz,deru(kr),derv(kr),derw(kr)
 92      format(1x,i5,e10.3,3e15.4)
      end do 
c     
 
 25   continue

c     
      do k=1,nz
         j=nyh
         i=nxh
         deru(k)= ox(i,j,k)
         derv(k)= oy(i,j,k)
         derw(k)= oz(i,j,k)
      end do
c     
      print *,'across vertical line: dudz, dvdz, dwdz ' 
       do k=1,nz
      kr=nz-k+1       
C     reverse index for printing from top to bottom
      zz=0.5*(zpts(kr-1)+1.0)
      write(6,92) kr,zz,deru(kr),derv(kr),derw(kr)
       enddo
                       
c    
c  go back to spectral space 

 35   continue
      call horfft(u,-1,nz)
      call horfft(v,-1,nz)
      call horfft(w,-1,nz)   
c
      return
       end         

c----------------------------------------------------------------------

C****************************************************************************
      subroutine tempreset(zlen)
C-Sets the temperature to zero for the beginning of a restart run.
C****************************************************************************
      include 'dim.h'

      include 'comflow.h'
      include 'comscratch.h'
      include 'comleng.h'
      include 'combc.h'
      include 'comgridptxyz.h'

       write (*,*) 'Resetting TEMPERATURE field upon restart.'
       call  horfft (temp,1,nz)
       call  horfft (ox,1,nz)
c-Set temperature values
       do k=1,nz
         zz = z(k)
         do j=1,ny
           do i=1,nxpp

           
c            if ((i.eq.nxh).and.(j.eq.nyh+1))
c    >       write(*,*) z(k),temp(i,j,k),tn(i,j,k),tnm(i,j,k)
             temp(i,j,k)=0. ! zz/zlen
             tn(i,j,k) = 0.
             tnm(i,j,k) = 0.
             ox(i,j,k) = 0.

c            if ((i.eq.nxh).and.(j.eq.nyh+1))
c    >       write(*,*) z(k),temp(i,j,k),tn(i,j,k),tnm(i,j,k)
           enddo
         enddo
       enddo
C-Reset boundary conditions
        kk = (nz-1)*nxppy
        do j=1,ny
          do i=1,nx
c-Fixed temperatures at boundaries
            ox(i,j,nz) = 0.0
            temp(i,j,nz)=0.0
            ox(i,j,1)  = 0.0
            temp(i,j,1)= 0.0
          enddo
        enddo

        call  horfft (temp,-1,nz)
        call  horfft (ox,-1,nz)
c
c     save the boundary condition along the upper plate (tubc)
c     and along the bottom plate (tbbc)
c
         do j=1,ny
           do i=1,nxpp
             tubc(i,j) = ox(i,j,nz)
             tbbc(i,j) = ox(i,j,1)
           enddo
         enddo
       
         return 
         end 

      





      subroutine hiter1(z2,nz2,z1,nz1,hf)
      dimension z2(nz2),z1(nz1),hf(nz2,nz1)
c     calculates the interpolation coeficients
c     z2 = mesh to be interpolated
c     z1 = original mesh
c     hf - array containing interpolation coeficient
       
      print *,'in sub, nz1,nz2',nz1,' --',nz2
      do k = 1,nz2
         do j = 1,nz1
            hf(k,j) = 1.0000
            do i = 1,nz1
               if (i.ne.j) then
                  hf(k,j) = hf(k,j)*(z2(k)-z1(i))/(z1(j)-z1(i))
               end if
            end do
         end do
      end do
      
      return
      end

c----------------------------------------------------      
      
      subroutine testdump(t,ybar,iflag,tstart)

C-a) Dumps out a sample profile at the center of the horizontal
C-plane.
C-b) Generates contour plots of u,v,w,T at center Oxz & Oyz planes.
C-(Only if iflag=1)
!
!-Set indices that determine rough boundaries of separation bubble to
!-perform averaging within.
      parameter (il=461,ir=691,kl=1,kr=107)

      include 'dim.h'
      include 'comflow.h'
      include 'commatrix.h'
      include 'comleng.h'
      include 'comwave.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'
      include 'comparamphys.h'
   
c************************************
      dimension xdudz(nxpp,ny,nz),xdwdx(nxpp,ny,nz),xdtdz(nxpp,ny,nz)
 
      common /deriv/ xdudz,xdwdx,xdtdz


      character*30 fileout,fileoutOxz,fileoutOyz,fileout2,testhobands
     >,fileout3,tmflow
      character*9 fmt1
      real, dimension(nz) ::  Ubar,Ufl,Ucor,Ubuoy
      real :: sumKE,sumPE,sumTE,sumbu,sumco,zz,Dengrad,
     > sum_mean,outpe,outke,outte,outhe,sumzpe,sumzte,sumzke,sumzhe
!    >        sumshear,Rhoder,Uzder,Dengrad 
      real, dimension(1:nz) :: XKE,XPE,XTE,XHE
      real, dimension(1:nzloc) ::ake,ape,ate,ahe
!-Set to maximum absolute and rms-averaged value of w-velocity in
!-separation bubble.
      fileout='testprof'
      fileoutOxz='cntroxz_'
      fileoutOyz='cntroyz_'
      fmt1 = '(f5.3,a4)'
      if (int(t).ge.10) fmt1 = '(f6.3,a4)'
      if (nint(t).ge.100) fmt1 = '(f7.3,a4)'
      if (nint(t).ge.1000) fmt1 = '(f8.3,a4)'
	if (nint(t).ge.10000) fmt1 = '(f9.3,a4)'
	
      write(unit=fileout(9:22),fmt=fmt1) t-tstart,'.dat'
      if (iflag.eq.1) then
        write(unit=fileoutOxz(9:18),fmt=fmt1) t-tstart,'.dat'
        write(unit=fileoutOyz(9:18),fmt=fmt1) t-tstart,'.dat'
      endif
      write(*,*) fileout 
     
      open(400,file=fileout)
      write(400,*) 
     >'VARIABLES = 
     >"z","u","v","w", "T","Ubar","Ufl","Ucor","Ubuoy"' 
     
      open(66,file='growthrate.dat',access='append') 

      if (iflag.eq.1) then 
        open(401,file=fileoutOxz,status='unknown')
        open(402,file=fileoutOyz,status='unknown')
        

        write(401,*) 'VARIABLES = "x","z","u","v","w","T"'
        write(402,*) 'VARIABLES = "y","z","u","v","w","T"'

        write(401,*) 'ZONE F=POINT, I=',nx,', J=',nz
        write(402,*) 'ZONE F=POINT, I=',ny,', J=',nz
      endif
       
      call horfft(u,1,nz)
      call horfft(v,1,nz)
      call horfft(w,1,nz)
      call horfft(temp,1,nz)

C-Create profile
C-(At horizontal position = dependent on type of wave forcing
!-and Reynolds number -- See timeseries generation below)

         open(698,file='xtu.dat',access='append')
         open(697,file='xtw.dat',access='append')
         open(696,file='xttemp.dat',access='append')

         open(700,file='uveltest.dat',access='append')
c       zmeasure is the approimate height at which horizontal
c       xt measurements are taken
        zmeasure = 0.175
        zmeasure2 = 0.2625
        zmeasure3 = 0.2433

        counter = 1
        switch = 0

        counter2 = 1
        switch2 = 0

        counter3 = 1
        switch3 = 0

        do while(switch.eq.0) 
           zcheck = z(counter)
           if(zcheck.lt.zmeasure) then
           counter = counter + 1
           else
           switch = 1
           endif
        enddo 

        do while(switch2.eq.0)
           zcheck2 = z(counter2)
           if(zcheck2.lt.zmeasure2) then
           counter2 = counter2 + 1
           else
           switch2 = 1
           endif
        enddo

        do while(switch3.eq.0)
           zcheck3 = z(counter3)
           if(zcheck3.lt.zmeasure3) then
           counter3 = counter3 + 1
           else
           switch3 = 1
           endif
        enddo

c        write(*,*) switch, counter, z(counter)      

        k = counter  ! int(9.5*float(nx)/xlen)
        k2 = counter2
        k3 = counter3
        if (ndims.eq.3) then 
          j = nyh
        else
          j = 1
        endif
!------------------------------------------------------
!-----------------------------------------------------      
        Dengrad= rhograd(zlen,zlen,grav,rho0,brunt)
!-----------------------------------------------------
        do i=1,nx

          delNsq = (-grav/rho0)*xdtdz(i,j,k)
          delNsq2 = (-grav/rho0)*xdtdz(i,j,k2)
          delNsq3 = (-grav/rho0)*xdtdz(i,j,k3)

         write(698,*) x(i),t,u(i,j,k),u(i,j,k2),u(i,j,k3)
         write(697,*) x(i),t,w(i,j,k),w(i,j,k2),w(i,j,k3)
         write(696,*) x(i),t,delNsq,delNsq2,delNsq3

        enddo
!---------------- Horizontal Averaging -------------------------
       open(738,file='ENERGYT.dat',access='append')
!          do k=1,nz 
!             XKE(k)=0.5*(horav(u,k,2)+horav(w,k,2))
!             XPE(k)=0.5*(horav(brunt*temp/Dengrad,k,2))
!             XHE(k)=0.5*(Ubar(k)+Ufl(k))**2 
!       XTE(k)=0.5*(Ufl(k)**2+2.*Ufl(k)*horav(u,k,1)+horav(u,k,2)
!     >    +horav(w,k,2)) 

!        end do
!--------------------------------
            XKE(1)=0.0; XKE(nz)=0.0
            XPE(1)=0.0;XPE(nz)=0.0
            XTE(1)=0.0;XTE(nz)=0.0
            XHE(1)=0.0;XHE(nz)=0.0            
!----------------------------------------------
       sumzke=0.0
       sumzpe=0.0
       sumzte=0.0
       sumzhe=0.0
      do ks=1,nsubd
        ybar = zh(ks)
        do kloc = 1,nzloc
          k = (ks-1)*nzloc + kloc
          ake(kloc)=XKE(k)
          ape(kloc)=XPE(k)
          ate(kloc)=XTE(k)
          ahe(kloc)=XHE(k)
        enddo
        call integrcol(ake,outke,ndim,ybar)
        sumzke = sumzke + outke
        call integrcol(ape,outpe,ndim,ybar)
        sumzpe=sumzpe + outpe
        call integrcol(ate,outte,ndim,ybar)
        sumzte=sumzte + outte
        call integrcol(ahe,outhe,ndim,ybar)
        sumzhe=sumzhe + outhe
        enddo
        write(738,123)t,sumzke/zlen,sumzpe/zlen,sumzhe/zlen,
     >  sumzte/zlen,(sumzte+sumzpe)/zlen  
C*************NOW AVERAGE IN THE Z************
!       write(738,123)t,verav(XKE),verav(XPE),verav(XHE),verav(XTE),
!     >     (verav(XTE)+verav(XPE))      
123    format (1x,f12.6,2x,5(e12.6,2x))
       close(738)
!---------------------------------------------------------
!---------------------------------------------------------
C-Create contour plot files
      if (iflag.eq.1) then

        do k=1,nz
          do i=1,nx
            write(401,'(2x,2(F10.6,1x),4(G15.8,1x))') x(i),z(k),
     > u(i,nyh,k),v(i,nyh,k),w(i,nyh,k),temp(i,nyh,k)
          enddo

          do j=1,ny
            write(402,'(2x,2(F10.6,1x),4(G15.8,1x))') y(j),z(k),
     >u(nxh,j,k),v(nxh,j,k),w(nxh,j,k),temp(nxh,j,k)
          enddo
        enddo

        close(401)
        close(402)

      endif

C-Dump out limit cycle info
!-Depression NLIW
!-Re=20K: at (k=51) z=0.175 and x=13.5, 14 & 14.5
!-Re=100K:  at (k =108 ) z=0.100 and x= 14, 14.5 & 15

!-Elevation NLIW
!-Re=100K: (at k=40) z=0.075 and x = 9, 9.5, & 10
!-  open(64,file='vertimeseries.dat')

        open(64,file='vertimeseries.dat',access='append')


       i1 =20    ! int(9.0*float(nx)/xlen)
       i2 =30 !int(9.5*float(nx)/xlen)
       i3 =32 !int(10.0*float(nx)/xlen)
       k0 =242
       k1=340

       write(64,'(1x,f12.6,2x,6(e12.6,2x))') t,u(i3,1,k0),w(i3,1,k0),
     >u(i3,1,k1),w(i3,1,k1)   !,u(i3,1,k0),w(i3,1,k0)
       close(64)
!
      write(66,'(1x,f12.6,2x,3(e12.6,2x))') t,wrms,wmax,zmax
          

      write(*,*)'malaka', t,u(32,1,10),u(32,1,58)
      call horfft(u,-1,nz)
      call horfft(v,-1,nz)
      call horfft(w,-1,nz)
      call horfft(temp,-1,nz)

      close(400)
      close(66) 
      
 100  format (1x,f12.6,2x,7(e12.6,2x))
    
      return
      end





C****************************************************************
      subroutine maxcalc(istep)
C****************************************************************
C-Calculates the minimum and maximum values of (u,v,w,T) and prints
C-out their k-indexed  location
      include 'dim.h'

      include 'comflow.h'
      include 'comgridptxyz.h'
      include 'comtime.h' 
      include 'comwave.h'
      include 'comsubd.h'
      include 'comsoliton.h'
      include 'comcritlayer.h'
      include 'comparamphys.h'
      include 'comtimeadvance.h'
     
c      include 'equivsuvwt.h'
c      include 'equivscratch.h'
c      include 'equivflow.h'
c     include 'equivflowprev.h'
c      include 'equivgrad.h'

 
      common /Dump/ Tdump,Tx

!      parameter (cflzmin_lim = 0.7, cflzmax_lim = 0.95,
!    >                   cflxmin_lim = 0.2)
      
      logical change_dt

C-Calculate coordinate values of grid-points and
C-store in corresponding common block
        pi2 = 8.0*atan(1.0)
        xl = pi2/alpha
        yl = pi2/beta
        dy = yl/ny
        dx = xl/nx      

        umin=0.
        umaxt=0.
        vmin=0.
        vmax=0.
        wmin=0.
        wmax=0.
        tmin=1.
        tmax=0.

        cflxmax = 0.
        cflymax = 0.
        cflzmax = 0.
 

        kumin=0
        kumax=0
        kvmin=0
        kvmax=0
        kwmin=0
        kwmax=0
        ktmin=0
        ktmax=0

        itmax=0
        jtmax=0

        icflx = 0
        icfly = 0 
        icflz = 0
        jcflx = 0
        jcfly = 0
        jcflz = 0
        kcflx = 0
        kcfly = 0
        kcflz = 0

        call  horfft (u,1,nz)
        call  horfft (v,1,nz)
        call  horfft (w,1,nz)
        call  horfft (temp,1,nz)

        open(700,file='uveltest.dat',access='append')
        write(700,*) t,u(32,1,10),u(32,1,50)
        close(700)

        do i=1,nx
          do j=1,ny
c           do k=1,nz

            do ks=1,nsubd
              do kloc = 1,nzloc
                k = (ks-1)*nzloc + kloc
                kg = (ks-1)*(nzloc-1) + kloc

            if ((u(i,j,k)+umwv(i,k)).ge.umaxt) then
              umaxt=u(i,j,k)+umwv(i,k)+umax
              kumax=k
            endif
            if ((u(i,j,k)+umwv(i,k)).le.umin)  then
              umin=u(i,j,k)+umwv(i,k)
              kumin=k
            endif

            if (v(i,j,k).ge.vmax) then
              vmax=v(i,j,k)
              kvmax=k
            endif
            if (v(i,j,k).le.vmin)  then
              vmin=v(i,j,k)
              kvmin=k
            endif

            if ((w(i,j,k)+wmwv(i,k)).ge.wmax) then
              wmax=w(i,j,k)+wmwv(i,k)
              kwmax=k
            endif
            if ((w(i,j,k)+wmwv(i,k)).le.wmin)  then
              wmin=w(i,j,k)+wmwv(i,k)
              kwmin=k
            endif

            zz = z(k)
            if ( (temp(i,j,k)+rhomwv(i,k)+
     >           rhomean(zz,zlen,grav,rho0,brunt)) .ge.tmax) then
              tmax= temp(i,j,k)+rhomwv(i,k)+
     >rhomean(zz,zlen,grav,rho0,brunt)
              ktmax=k
              itmax=i
              jtmax=j
            endif

            if ( (temp(i,j,k)+rhomwv(i,k)+
     >           rhomean(zz,zlen,grav,rho0,brunt)) .le.tmin)  then
              tmin= temp(i,j,k)+rhomwv(i,k)+
     >rhomean(zz,zlen,grav,rho0,brunt)
              ktmin=k
              itmin=i
              jtmin=j
            endif

!           advmean = umwv(i,k) - umax +
c              if(t.lt.600)fcshear=0
c              if(t.ge.0)fcshear=1
c             fcshear=1.0             
              
             cflx=abs((uo*fcshear*fshear(zz,zdim,zsh,xacu)+u(i,j,k))
     >                 *dt/dx)
c              write (*,*)'uo=',uo,'fcshear=',fcshear,'cflx=',cflx
c              write(*,*)'zdim',zdim,'zsh',zsh,'xacu',xacu

!    >      abs(umax)*fcshear*fshear(zz,zlen) 
c            cflx = abs(u(i,j,k)*dt/dx)
            cfly = abs(v(i,j,k)*dt/dy)
C-The following few lines should normally be modified
C-to incorporate soliton velocities at top boundary.
            if ((k.eq.0).or.(k.eq.nz)) then
              cflz = 0.
            else
              dz = (zf(kg)-zf(kg-1))
              cflz = abs((w(i,j,k)*dt)/dz)
            endif

            if (cflx.ge.cflxmax) then
              cflxmax = cflx
              icflx = i
              jcflx = j
              kcflx = k
            endif

            if (cfly.ge.cflymax) then
              cflymax = cfly
              icfly = i
              jcfly = j
              kcfly = k
            endif

            if (cflz.ge.cflzmax) then
              cflzmax = cflz
              icflz = i
              jcflz = j
              kcflz = k 
              dzmax = dz
            endif

          enddo
          enddo
c         enddo
         enddo
        enddo


        call  horfft (u,-1,nz)
        call  horfft (v,-1,nz)
        call  horfft (w,-1,nz)
        call  horfft (temp,-1,nz)

        write(50,
     +'(1x,I5,2x,4(E12.6,1x,E12.6,1x,I3,1x,I3,2x))')
     +  istep,umin,umaxt,kumin,kumax,
     +        vmin,vmax,kvmin,kvmax,
     +        wmin,wmax,kwmin,kwmax,
     +        tmin,tmax,ktmin,ktmax
c    +'(1x,I5,2x,e12.6,2x,e12.6,2x,i2,2x,i2,1x)')
c    + istep,tmin,tmax,kmin,kmax
        write(*,
     +'(1x,I5,2x,e12.6,2x,e12.6,2x,i3,2x,i3)')
     + istep-1,tmin,tmax,ktmin,ktmax

        write(51,'(1x,e12.6,2x,3(1x,e12.6,2x,3(I5,2x)))')
c       write(*,*)
     >t,
     >cflxmax,icflx,jcflx,kcflx,
     >cflymax,icfly,jcfly,kcfly,
     >cflzmax,icflz,jcflz,kcflz


       write(*,*) 'Temp. min at i,j,k',itmin,jtmin,ktmin
       write(*,*) 'Temp. max at i,j,k',itmax,jtmax,ktmax

       write(*,*) 'At max CFLZ dz= ',dzmax, cflzmax
!-Adapt Timestep if CFL# in z-direction lies outside
!-allowed range of [0.25,1] (see beginning of this routine)
!-AND (CAREFUL) if timestep does not exceeded the maximum
!-timestep imposed by the stratification period.
!-cdt indicates how much bigger older timestep is w/r
!-to new one.
!-Note: Adaptive timestepping scheme applies only if we have
!-crossed t > 1 sec.
c       write(*,*)'time=',t,'xomega=',xomega,'Period=',pi2/xomega
        cflxmax_lim=0.16
 
c           cflzmax_lim=1.e8
c           cflzmin_lim=1.e-8


c         if (t.lt.pi2/xomega) then
c            cflzmax_lim=1.e9        
c            cflzmin_lim=1.e-9
c          else 
            cflzmax_lim=0.7
            cflzmin_lim=0.1
c          endif

c       change_dt = ( (cflzmax .gt. cflzmax_lim).or.
c     >      ((cflzmax .lt. cflzmin_lim).and.
c     >       (cflxmax .lt. cflxmin_lim))
c     >        .and.(int(t-tstart).gt.0) )
c%************************************************************
           change_dt = ( (cflzmax .gt. cflzmax_lim).or.
     >    (cflzmax .lt. cflzmin_lim).or.(cflxmax .gt. cflxmax_lim))

              
C****************ADAPTIVE TIME STEPPING. A.M.A. 06-11-08******
!-PD-Lyon-7/4/15: Changed Ammar's approach
         Tlim=pi2/(xomega*50.) ! Max. time step imposed by WAVE FREQUENCY (PD)
         Nst=(t-Tso)/dt   ! number of time steps since last time step change

         if (Nst.gt.20)then 

C**********************************************************
         if (change_dt) then
            Tso=t
c**********************************************************           
C********** The max function to  account for a case in which
C*********** dt is smaller than Tlim however when it is increased
C********** by a factor of 5/4 it becomes bigger than Tlim
C********** As when dt/Tlim > 0.8 this means dt/0.8 > Tlim 
C********** and we do not allow that.

          if (cflzmax .lt. cflzmin_lim) cdt=max(0.8,dt/Tlim)
c*************************************************************** 
          if ((cflzmax .gt. cflzmax_lim).or.
     >       (cflxmax.gt.cflxmax_lim)) cdt=1.25     
C******************************************************************
c         if (cflzmax .gt.cflzmax_lim) cdt=1.25
 
c            write(*,*)'cdt=',cdt

          dt = dt/cdt

C-Set up dt1,dt2,dt3 (and time-advancement coefficients).
C-Rationale same as restart. See main.f subroutine.
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

C-Adjust facvel, fact
         fact = 1./(xkappa*dt)
         facvel = 1./(pr*dt)

C-Set flags associated with adaptive timestepping
         ivardt = 1
         idtstep = 0
          write(65,*)t,dt
c         write(65,'(A10,F8.3,A21,F8.3)') 'At time t=',t,
c     > ', timestep set to Dt=',dt
c        write(65,'(A14,F8.3)') 'Because CFLZ=',cflzmax
c         write(65,'(3(F8.3))')t,dt,cflzmax
         write(*,*) 'At time= ',t
         if (cdt .gt. 1.)
     >   write(*,*) 'TIMESTEP *REDUCED* by a factor of 4/5 to dt=',dt
         if (cdt .lt. 1.)
     >   write(*,*) 'TIMESTEP *INCREASED* by a factor of 5/4 to dt=',dt
c         if (cdt .eq. 1.)
c     >   write(*,*) 'TIMESTEP *FIXED* at Max Set by Stratification 
c     >               dt=',dt
c         write(65,*) 'Factor cdt=',cdt
      
         endif

        endif 

       return
       end

C-LINUX Timer routines (Installed by PD 2/22/05)

******************************************************************
       subroutine timecall(isecout)
C******************************************************************
       character*9 ctime
       integer isec,imin,ihour
       ifac = 48
                                                                                
       call time(ctime)
       ihour = (ichar(ctime(1:1)) - ifac)*10. +
     >        ichar(ctime(2:2)) - ifac
       imin = (ichar(ctime(4:4)) - ifac)*10. +
     >        ichar(ctime(5:5)) - ifac
       isec = (ichar(ctime(7:7)) - ifac)*10. +
     >        ichar(ctime(8:8)) - ifac
                                                                                
c-Now compute total seconds elapsed
       isecout = ihour*3600+
     >           imin*60 + float(isec)
                                                                                
       return
       end

C******************************************************************
       subroutine timeelapse(isec0,isec1)
C******************************************************************
       timetot = float(isec1-isec0)
       write(*,*) 'Time Elapsed for this operation in seconds',timetot
     >
       return
       end

C-Timer routines that exploit intrinsic wrappers in UNIX.
C-(PD: 2/22/05) Do not use them when running on Linux.
C******************************************************************
!      subroutine timecall(isec)
C******************************************************************
!      integer isec,time
!      isec = time()
!      return
!      end

C******************************************************************
!      subroutine timeelapse(isec0,isec1)
C******************************************************************
!      timetot = float(isec1-isec0)
!      write(*,*) 'Time Elapsed for this operation in seconds',timetot
!      return
!      end

C-Subroutines below employed IMSL time calls. IMSL was jinxed
C-not to interfere with BLAS (PD: 7/25/03) 
C******************************************************************
c      subroutine timecall_old(ihour,imin,isec)
C******************************************************************
C-Calculates time in hrs,mins,secs via imsl routine.
c      integer ihour,imin,isec,nout
c      external timdy,umach

c      call timdy(ihour,imin,isec)
c      call umach(2,nout)

c      return
c      end

C******************************************************************
c      subroutine timeelapse_old(ihour0,ihour1,imin0,imin1,isec0,isec1)
C******************************************************************
C-Calculates time elapsed between two consecutive calls of timdy.
c      integer ihour0,imin0,isec0,ihour1,imin1,isec1

c      idh = ihour1-ihour0
c      idm = imin1-imin0
c      ids = isec1-isec0

c      timetot = float(idh)*3600. + float(idm)*60. + ids

c      write(*,*) 'Time Elapsed for this operation in seconds',timetot

c      return
c      end

       

