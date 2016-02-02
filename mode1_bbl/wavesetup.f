c-This file contains all necessary operations
C-that set up the forcing of a solitary wave
C-of given: wavenumber, mode

c*********************************************************************
      subroutine wavesetup
c*********************************************************************
      include 'dim.h'

      include 'comflow.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comsoliton.h' 
      include 'comcritlayer.h'

      dimension f(nzloc),df(nzloc),d2f(nzloc)
      dimension fin(nzloc),fout(nzloc)
      dimension testf(nz)

C-Calculate dimensional amplitude of A(x,t) function
C-Coefficient alpha0
C-Non-Sheared Case
       alpha0=(zeta0/zlen)*abs(umax)*zlen
C-signs of alpha0
      if (alphak.gt.0) then
        a0 = -alpha0
        write(*,*) 'Wave of Elevation'
      else
        a0 = alpha0
        write(*,*) 'Wave of Depression'
      endif

       pi2 = 8.*atan(1.0)
       xl = pi2/alpha
       dx = xl/nx

      do ks=1,nsubd
        do kloc=1,nzloc
          k = (ks-1)*nzloc + kloc
          kg = (ks-1)*(nzloc-1) + kloc
          dz = (zf(kg)-zf(kg-1))
          zz = z(k)
          do i=1,nx
            xx = x(i)
C-Velocity
            umwv(i,k) = 0.0
            wmwv(i,k) = 0.0
C-Velocity Gradients (1st Derivatives)
            dumdx(i,k) = 0.0
            dumdz(i,k) = 0.0                 
            dwmdx(i,k) = 0.0
            dwmdz(i,k) = 0.0
C-Density   
          

            c1 = (umax*(1. - fcshear*fshear(z,zdim,zsh,xacu)
     >           *float(isignc)))
 
            rhomwv(i,k) = 0.0
            drhomdx(i,k)= 0.0
            drhomdz(i,k)= 0.0

C-Now dump data to file
            drhomean = abs( rhomean(z(nz),zlen,grav,rho0,brunt) -
     >               rhomean(z(1),zlen,grav,rho0,brunt) )

            write(600,'(2x,6(e12.6,1x))')
     >      xx,zz,umwv(i,k),wmwv(i,k),rhomwv(i,k) 
          enddo
        enddo
      enddo

c     close(600)

C-Dump out sample soliton velocity profile at middle
C-of domain
      k = nz
c     i = 240
      open(601,file='soliton.dat')
      write(601,*) 'VARIABLES = "x","`z","u","w","dudz","dwdz"'
      do i=1,nx
c     do k=1,nz
        xx=x(i)
        zz=z(k)

c         write(601,'(2x,6(e12.6,1x))') zz,fsol(zz,zlen,nmode),
c    >    dfsoldz(zz,zlen,nmode),
c    >    d2fsoldz2(zz,zlen,nmode),
c    >    d3fsoldz3(zz,zlen,nmode),
c    >    fsol(zz,zlen,nmode)
 
          write(601,'(2x,6(e12.6,1x))') xx,asol(xx,xlen,solks,a0),
     >    dasoldx(xx,xlen,solks,a0),
     >    d2asoldx2(xx,xlen,solks,a0),
     >    d3asoldx3(xx,xlen,solks,a0),
     >    asol(xx,xlen,solks,a0)

c    >    umwv(i,k),wmwv(i,k),dumdx(i,k),dwmdx(i,k)
c         write(601,'(2x,6(e12.6,1x))') xx,zetamwv(xx,xlen,zeta0,solks),
c    >    umwv(i,k),wmwv(i,k),dumdx(i,k),dwmdx(i,k)

      enddo
      close(601)

      return
      end
