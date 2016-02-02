c**********************************************************************
      subroutine output(nhist,t)                                         
c**********************************************************************
c                                                                       
      include 'dim.h'
c    
c      
      include 'comwave.h'
      include 'comleng.h'
      include 'comflow.h'
      include 'comcoeff.h'
      include 'comparamphys.h'
      include 'combc.h'
      include 'comlegendre.h'
      include 'comfilter.h'
      include 'comddzwfunc.h'
      include 'comzeta.h'
      include 'comwfunc.h'
      include 'comsubd.h'
      include 'comgridptxyz.h'
      include 'commatrix.h'
      include 'comsoliton.h'
            
      
      common /Dump/ Tdump,Tx
      common /anim/ianim,ianimc 
C-Dumps field at existing and previous timestep. Also a bunch
C-of auxilliary information generated in subroutine setup &
C-elsewhere. Thus, one does not need to run setup upon restart
      write(nhist) t 
      write(nhist) pr,ra
      write(nhist) nx,ny,nz,nzm,nxppy,ntot,nsubd,nzloc,nsdtype
      write(nhist) alpha,beta                                            
      write(nhist) (xw(i),xsq(i),i=1,nx)                             
      write(nhist) (yw(i),ysq(i),i=1,ny)                             
      write(nhist) (x(i),i=1,nx)
      write(nhist) (y(j),j=1,ny)
      write(nhist) (z(k),k=1,nz)
      write(nhist) (zf(k),k=1,nz2)
      write(nhist) (z0(k),zh(k),isubdtype(k),k=1,nsubd)
      write(nhist) (subdlen(k),k=1,nsdtype)
      write(nhist) ((alpol(j,k),j=0,nzm),k=0,nzm)
      write(nhist) ((dalpol(j,k),j=0,nzm),k=0,nzm)
      write(nhist) (zpts(i),i=0,nzm)                                     
      write(nhist) (wg(i),i=0,nzm)
      write(nhist) ((d(i,j),i=0,nzm),j=0,nzm)
      write(nhist) ((d2(i,j),i=0,nzm),j=0,nzm)
      write(nhist) ((d3(i,j),i=0,nzm),j=0,nzm)
      write(nhist) (atpz(k),btpz(k),ctpz(k),k=1,nz)
C************************************************************************
c      write(nhist) ((tubc(i,j),tbbc(i,j),i=1,nxpp),j=1,ny)
       write(nhist) ((tubc(i,j),i=1,nxpp),j=1,ny)
       write(nhist) ((tbbc(i,j),i=1,nxpp),j=1,ny)
c      write(nhist) (((uw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2),
c     &             (((vw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2),
c     &             (((ww(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2)
C*** I read vars separately, more convenient for regridding purposes
C** A.M.A. Oct-13-08
        write(nhist) (((uw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2)
        write(nhist) (((vw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2) 
        write(nhist) (((ww(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2)

c Soliton velocity/density and gradient fields (includes pressure gradient)
c      write(nhist) solks
c      write(nhist) ((umwv(i,k),wmwv(i,k),i=1,nxpp),k=1,nz)
c      write(nhist) ((dumdx(i,k),dumdz(i,k),i=1,nxpp),k=1,nz)
c      write(nhist) ((dwmdx(i,k),dwmdz(i,k),i=1,nxpp),k=1,nz)
c      write(nhist) ((rhomwv(i,k),i=1,nxpp),k=1,nz)
c      write(nhist) ((drhomdx(i,k),drhomdz(i,k),i=1,nxpp),k=1,nz)
       
      
c velocity and temperature at (n) time level in spectral-physical space
c
      write(nhist) ( ((u(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((v(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((w(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((temp(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)

c
c velocity and temperature (n-1) time level in physical space
c (for BDF purposes)
c
      write(nhist) ( ((un(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((vn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((wn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((tn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
c
c
c velocity and temperature (n-2) time level in physical space
c (for BDF purposes)
c
      write(nhist) ( ((unm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((vnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((wnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((tnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
 

c
c convective terms at (n-1) time level in physical space
c
      write(nhist) ( ((sun(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((svn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((swn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((stn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
c
c
c convective terms at (n-1) time level in physical space
c
      write(nhist) ( ((sunm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((svnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((swnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      write(nhist) ( ((stnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
c 
      write(nhist) ianimc,Tdump,Tx
      write(6,*) 'finished writing to tape ',nhist 
      write(6,*) 'at time ',t
      
c     do k =1,nz
c       i = nxh  
c       j = nyh+1
c       write(*,'(5(E12.6,2x))') 
c    >  z(k),u(i,j,k),v(i,j,k),w(i,j,k),temp(i,j,k)
c     enddo
c
      return                                                            
      end                                                               






c**********************************************************************
      subroutine plotf(iunit,iform,t)   
c**********************************************************************
c                                                                       
      include 'dim.h'
cc
      include 'comwave.h'
      include 'comleng.h'
      include 'comflow.h'
      include 'comgridptxyz.h'

C-PD: 7/28/03.
C-The purpose of tempout (...) is to serve as a temporary
C-4-byte real array to which output data of (u,v,w,T) is
C-temporarily copied over before output. Thus, one consumes
C-much less disk-space for 3-D field output. Make sure that
C-the separate postprocessing program is aware of this.
      real*4 tempout(nxpp,ny,nz)
             
c     dimension  ud(*),vd(*),wd(*),td(*)                 
c          ,u(ntot),
c    &     v(ntot),w(ntot),temp(ntot)
      
c
c
C-Convert to Physical Space for output
      call horfft(u,1,nz)
      call horfft(v,1,nz)
      call horfft(w,1,nz)
      call horfft(temp,1,nz)
c
      

      if(iform.eq.0) then
         write(iunit) t            
         write(iunit) nx,ny,nz,ntot
         write(iunit) (x(i),i=1,nx)
         write(iunit) (y(i),i=1,ny)
         write(iunit) (z(i),i=1,nz)
C-The temporary copies are not done through some subroutine
C-to avoid screw-ups with memory addressing (PD: Improve this asap !)
C-Dump u-velocity 
         do k=1,nz
           do j=1,ny
             do i=1,nxpp
 
               tempout(i,j,k) = real(u(i,j,k))
 
             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)

C-Dump v-velocity 
         do k=1,nz
           do j=1,ny
             do i=1,nxpp
  
               tempout(i,j,k) = real(v(i,j,k))

             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
        
C-Dump w-velocity 
         do k=1,nz
           do j=1,ny
             do i=1,nxpp
  
               tempout(i,j,k) = real(w(i,j,k))

             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,ny),k=1,nz) 

C-Dump temperature
         do k=1,nz
           do j=1,ny
             do i=1,nxpp
  
               tempout(i,j,k) = real(temp(i,j,k))

             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
c     
      endif
      
      if(iform.eq.1) then 
         write(iunit,100) t            
         write(iunit,110) nx,ny,nz,ntot
         write(iunit,100) (x(i),i=1,nx)
         write(iunit,100) (y(i),i=1,ny)
         write(iunit,100) (z(i),i=1,nz)
         write(iunit,100) ( ((u(i,j,k),i=1,nxpp),j=1,ny),k=1,nz),
     >                ( ((v(i,j,k),i=1,nxpp),j=1,ny),k=1,nz),
     >                ( ((w(i,j,k),i=1,nxpp),j=1,ny),k=1,nz),
     >                ( ((temp(i,j,k),i=1,nxpp),j=1,ny),k=1,nz) 

      endif
c     

c //// to dump it to matlab file 
c
c  iform=10
c
c ////////////////////////////////
     
c      if(iform.eq.10) then 
c         do k=1,nz
c            do j=1,ny   
c               do i=1,nx
c       write(iunit,1003) i,j,k,u(i,j,k),v(i,j,k),w(i,j,k),temp(i,j,k),
c    +                 x(i),y(j),z(k)      
c               enddo
c            enddo 
c         enddo
c 1003    format(1x,3(i5,1x),1x,7(e9.3,1x),1x)
c         
c      endif
c ////////////////////     

C-PD: Reconvert to Fourier space now that output has been completed
      call horfft(u,-1,nz)
      call horfft(v,-1,nz)
      call horfft(w,-1,nz)
      call horfft(temp,-1,nz)
      

 100  format(5e14.6)
 110  format(5i10)
c
      return
      end

