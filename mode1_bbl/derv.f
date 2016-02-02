C-All routines in this file either calculate divergence or vertical
C-derivatives

c*********************************************************************  
      subroutine div(u,v,w,xw,yw,su)
c*********************************************************************  
c
c   calculate the divergence of a velocity field
c
      include 'dim.h'
c
      include 'comleng.h'
      include 'comscratch.h'
      dimension w(nxpp,ny,nz),xw(*),yw(*)      

      complex o3(nxhp,ny,nz)
      complex su(nxhp,ny,nz),u(nxhp,ny,nz),v(nxhp,ny,nz)
      equivalence (oz,o3)
c
      call dpdz (w,oz)
c
       do k=1,nz
        do j=1,ny
         do i=1,nxhp
           su(i,j,k) = (0.,1.)*(xw(i)*u(i,j,k) + yw(j)*v(i,j,k)) +
     1                  o3(i,j,k) 
         enddo
        enddo
       enddo  
       return
       end 
c

c*********************************************************************
      subroutine divcolc(ucol,vcol,wcol,xw,yw,divvel,dudx,dvdy,dwdz)
c*********************************************************************
c
c   calculate the divergence of a velocity field
c   for a column of data in Fourier-spectral space
c   (**PD: 1/3/03**)
c
      include 'dim.h'
c
      include 'comleng.h'
      include 'comsubd.h'
      complex rx
      complex ucol(nz),vcol(nz),wcol(nz),divvel(nz)
      complex dudx(nz),dvdy(nz),dwdz(nz)
      dimension fi(nzloc),fr(nzloc),dfi(nzloc),dfr(nzloc)
      rx = -1.0*cmplx(0.,1.)

c
C-Scan over all subdomains
       do 100 ks=1,nsubd
         ybar = zh(ks)
C-Assign values to temporary arrays
         do kloc = 1,nzloc
           k = (ks-1)*nzloc + kloc
           fr(kloc) = real(wcol(k))
           fi(kloc) = real(rx*wcol(k))
         enddo

C-Calculate vertical derivative in specific subdomain
         call dpdzcol(fr,dfr,d,zpts,nzm,ybar,1,nzloc) 
         call dpdzcol(fi,dfi,d,zpts,nzm,ybar,1,nzloc)

C-Now calculate divergence in specific subdomain
         do  kloc=1,nzloc
            k = (ks-1)*nzloc + kloc

C-Add d/dx and d/dy contributions
          dudx(k) = cmplx(0.,xw)*ucol(k)
          dvdy(k) = cmplx(0.,yw)*vcol(k)
          divvel(k) = divvel(k) + dudx(k) + dvdy(k)

C-Add d/dz contribution
          dwdz(k) = cmplx(dfr(kloc),dfi(kloc))
          divvel(k) = divvel(k) + dwdz(k)

c         write(*,*) dudx(k),dvdy(k),dwdz(k)

         enddo          
 100   continue
       return
       end
c


c*********************************************************************  
      subroutine divc(u,v,w,xw,yw,su)
c*********************************************************************  
c
c   calculate and check the divergence of a velocity field
c
      include 'dim.h'
c
      include 'comleng.h'
      include 'comscratch.h'

      dimension w(ntot),xw(*),yw(*)      

      complex o3(nxhp,ny,nz)
      complex su(nxhp,ny,nz),u(nxhp,ny,nz),v(nxhp,ny,nz)
      equivalence (oz,o3)
c
       call copy(w,o3,ntot)
       do i=1,nxh,4
       do j=1,ny,4
        write(6,12) i,j
 12     format(1x,'w at (i,j) location =',2i5)
       do k=1,nz
       write(6,13) k,o3(i,j,k)
 13    format(1x,i4,3x,2e12.3)
       end do
        end  do
         end do

c
c      call dpdz (w,oz,d,nzm,1.0)
c
       do k=1,nz
        do j=1,ny
         do i=1,nxhp
           u(i,j,k) = (0.,1.)*xw(i)*u(i,j,k)
           v(i,j,k) = (0.,1.)*yw(j)*v(i,j,k)
           su(i,j,k) = u(i,j,k) + v(i,j,k) + o3(i,j,k) 
         enddo
        enddo
       enddo
c
       do i=1,nxh,4
       do j=1,ny,4
        write(6,11) i,j
 11     format(1x,'(i,j) location =',2i5)
       do k=1,nz
       write(6,10) k,u(i,j,k),v(i,j,k),o3(i,j,k),su(i,j,k)
 10    format(1x,i4,4(3x,2e12.3))
       end do
        end  do
         end do
c  
       return
       end 


c********************************************************************
      subroutine dpdz_slow(f,df)
c*********************************************************************
C-Acts as an intermediary between the subroutine
C-that calls it and dpdz_slab which performs the actual
C-calculation.
      include 'dim.h'

      include 'comsubd.h'
      include 'comleng.h'

c     include 'comindex.h'
c
c     dimension a(ntot),b(ntot),d(0:ndim,0:ndim)
      dimension f(nxpp,ny,nz),df(nxpp,ny,nz)

      dimension temp(nxpp,ny,nzloc),dtemp(nxpp,ny,nzloc)
c
C-Scan over all subdomains
      do 100 ks = 1,nsubd

C-Transfer slab of (i,j,k) data associated with this subdomain
C-to temporary array.
        do kloc = 1,nzloc
          do j = 1,ny
            do i=1,nxpp
              
              k = (ks-1)*(nzloc-1) + kloc 
              temp(i,j,kloc) =  f(i,j,k)

c             if ((i.eq.1).and.(j.eq.1)) write(*,*) temp(i,j,kloc)
 
            enddo
          enddo
        enddo

C-Call routine for calculation of derivatives in slabs
C-NOTE: Given that we have C_1 continuity (Derivative continuity)
C-at the subdomain interfaces, the values of derivatives at
C-the bottom boundary of one subdomain should be the same
C-with that the top boundary of the subdomain below it.
         ybar = zh(ks)
C        call dpdz_slab(temp,dtemp,d,nzm,ybar)
C-Now re-update corresponding global derivative array
         do kloc = 1,nzloc 
          do j = 1,ny
            do i=1,nxpp
               
              k = (ks-1)*(nzloc-1) + kloc
              df(i,j,k) = dtemp(i,j,kloc) 
 
            enddo
          enddo
        enddo

c       write(*,*) 'Finished subdomain',ks

 100  continue

      return
      end
c



c********************************************************************
      subroutine dpdz(a,b)
c*********************************************************************
C-This is the old version of dpdz. A slab is taken out of the 3-D
C-array containing the variable of interest and the derivatives
C-are calculated within that slab (i.e. within the corresponding
C-spectral subdomain).
      include 'dim.h'
      include 'comsubd.h'
      include 'comleng.h'

c     include 'comindex.h'
c
c     dimension a(ntot),b(ntot),d(0:ndim,0:ndim)
      dimension a(nxpp,ny,nz),b(nxpp,ny,nz) ! d(0:nzm,0:nzm)
      dimension sum(nxpp,ny)
c
c finding the derivative of a function in the z-direction 
c
C-Scan over all subdomains
      do 110 ks=1,nsubd
C-Set domain height
        ybar = zh(ks)
        fac = 2./ybar

c
C-For a given local z-coordinate in subdomain, scan
C-over all x,y points.
        do k1=1,nzloc
          kk1 = (ks-1)*nzloc + k1
   
C-Nulify sum
          do i=1,nxpp
            do j=1,ny
              sum(i,j) = 0.0
            enddo
          enddo

C-Scan over all z-points in the subdomain to calculate
C-their contribution to the k1th point

          do k2=1,nzloc 
            kk2 = (ks-1)*nzloc + k2
c
            do i=1,nxpp
              do j=1,ny
                sum(i,j) = sum(i,j) + d(k1-1,k2-1)*a(i,j,kk2)
               enddo
            enddo
c
          enddo
c
          do i=1,nxpp
            do j=1,ny
              b(i,j,kk1) = sum(i,j)*fac
            enddo
          enddo
c
        enddo
c
  110 continue

c
      return
      end
c

c********************************************************************
      subroutine dpdzcol(a,b,d,zpts,ndim,ybar,nzbot,nztop)
c*********************************************************************
c-Using properties of Legendre interpolants calculates 1ST ORDER
c-vertical derivative of m-th order
c-within a specified range of datapoints in  a vertical column of data 
C-(6/20/02) WITHIN A GIVE VERTICAL SUBDOMAIN 

c
c-a: The input column of data
c-b: The output column. Only specified points are updated.
c-nzbot: Bottom point of point range
c-nztop: Top point of point range
      include 'dim.h'
      include 'comgridptxyz.h'

      dimension zpts(0:ndim)

c     include 'comindex.h'
c
c     dimension a(ntot),b(ntot),d(0:ndim,0:ndim)
      dimension a(*),b(*),d(0:ndim,0:ndim)
c     dimension a(nzloc),b(nzloc),d(0:ndim,0:ndim)
c
C-Set up values of collocation points

c
c-k1: (k-plane when scanning whole dataset, with the old indexing approach)
c-z-position we're interested in
      do k1=nzbot,nztop
c     do k1=nztop,nzbot,-1
 
          sum = 0.0
c
c-k2 index to scan z direction in
        do k2=1,nzloc
c       do k2=nzloc,1,-1
c
            sum = sum + d(k1-1,k2-1)*a(k2)
c
        enddo
c
c-Multiply by mapping factor. Z-location is taken into account
c-for the possibility of a more generalized mapping  
C-NOTE: In spectral subdomain decomposition, mapping not used !
          coeffmap = zderivmap(ybar)
          b(k1) = sum*coeffmap
c
      enddo
c
      return
      end
c




c********************************************************************
      subroutine d2pdz2col(a,b,d,zpts,ndim,ybar,nzbot,nztop)
c*********************************************************************
c-Using properties of Legendre interpolants calculates 1ST ORDER
c-vertical derivative of 2nd order
c-within a specified range of datapoints in  a vertical column of data
C-(6/20/02) WITHIN A GIVE VERTICAL SUBDOMAIN 
C-PD: 2/14/03. This could have been designed a little more efficiently.
C-I.e. a common routine for all differentiations. Just enter the order
C-m and the degree of differentiation. To avoid any mistakes, I stuck
C-with this.
c
c
c-a: The input column of data
c-b: The output column. Only specified points are updated.
c-nzbot: Bottom point of point range
c-nztop: Top point of point range
      include 'dim.h'
      include 'comgridptxyz.h'

      dimension zpts(0:ndim)

c     include 'comindex.h'
c
c     dimension a(ntot),b(ntot),d(0:ndim,0:ndim)
      dimension a(*),b(*),d(0:ndim,0:ndim)
c
C-Set up values of collocation points

c
c-k1: (k-plane when scanning whole dataset, with the old indexing approach)
c-z-position we're interested in
      do k1=nzbot,nztop

          sum = 0.0
c
c-k2 index to scan z direction in
        do k2=1,nzloc
c
            sum = sum + d(k1-1,k2-1)*a(k2)
c
        enddo
c
c-Multiply by mapping factor. Z-location is taken into account
c-for the possibility of a more generalized mapping
C-NOTE: In spectral subdomain decomposition, mapping not used !
          coeffmap = zderivmap2(ybar)
          b(k1) = sum*coeffmap
c
      enddo
c
      return
      end


c********************************************************************
      subroutine d3pdz3col(a,b,d,zpts,ndim,ybar,nzbot,nztop)
c*********************************************************************
c-Using properties of Legendre interpolants calculates 1ST ORDER
c-vertical derivative of 3rd order
c-within a specified range of datapoints in  a vertical column of data
C-(6/20/02) WITHIN A GIVE VERTICAL SUBDOMAIN
C-PD: 2/14/03. This could have been designed a little more efficiently.
C-I.e. a common routine for all differentiations. Just enter the order
C-m and the degree of differentiation. To avoid any mistakes, I stuck
C-with this.
c
c
c-a: The input column of data
c-b: The output column. Only specified points are updated.
c-nzbot: Bottom point of point range
c-nztop: Top point of point range
      include 'dim.h'
      include 'comgridptxyz.h'

      dimension zpts(0:ndim)

c     include 'comindex.h'
c
c     dimension a(ntot),b(ntot),d(0:ndim,0:ndim)
      dimension a(*),b(*),d(0:ndim,0:ndim)
c
C-Set up values of collocation points

c
c-k1: (k-plane when scanning whole dataset, with the old indexing approach)
c-z-position we're interested in
      do k1=nzbot,nztop

          sum = 0.0
c
c-k2 index to scan z direction in
        do k2=1,nzloc
c
            sum = sum + d(k1-1,k2-1)*a(k2)
c
        enddo
c
c-Multiply by mapping factor. Z-location is taken into account
c-for the possibility of a more generalized mapping
C-NOTE: In spectral subdomain decomposition, mapping not used !
          coeffmap = zderivmap3(ybar)
c         write(*,*) coeffmap,sum
          b(k1) = sum*coeffmap
c
      enddo
c
      return
      end

c*************************************************************************
      subroutine laplacian(finr,fini,flaplr,flapli,ak,ybar)
c*************************************************************************
C-Calculates the Laplacian of a complex function f(finr,fini) 
C-(vertical column of data) for
C-a given wavenumber ak in spectral space
      
       include 'dim.h'

       include 'comleng.h'
       include 'comcoeff.h'
    
       dimension finr(nzloc),fini(nzloc),fauxr(nzloc),fauxi(nzloc),
     + flaplr(nzloc),flapli(nzloc)

c-Nulify ? (No, not necessary)
       
c-first calculate 2nd derivative

       call d2pdz2col(finr,flaplr,d2,zpts,nzm,ybar,1,nzloc)

       call d2pdz2col(fini,flapli,d2,zpts,nzm,ybar,1,nzloc)

c-Update with contribution of horizontal wavenumbers
       do k=1,nzloc
      
         flaplr(k) = flaplr(k) - ak*finr(k)
         flapli(k) = flapli(k) - ak*fini(k)

       enddo

       return
       end
       
c********************************************************************
      subroutine integrcol(a,outint,ndim,ybar)
c*********************************************************************
c-Using properties of Legendre interpolants calculates 
c-integral of a given function in interval of length ybar.
C-ONLY WITHIN A GIVE VERTICAL SUBDOMAIN
C-(Created on 4/2/04 by PD)
c
c
c-a: The input column of data
c-b: The output value for integral
c-nzbot: Bottom point of point range
c-nztop: Top point of point range
      include 'dim.h'
      include 'comleng.h'
      include 'comgridptxyz.h'

      dimension a(*)
c
C-Set up values of collocation points

c
c-k: Local GLL grid-point
      sum = 0.0

      do k=1,nzloc
          sum = sum + wg(k-1)*a(k)
      enddo

c     write(*,*) sum,ybar
      outint = sum*ybar/2.
c
      return
      end             




