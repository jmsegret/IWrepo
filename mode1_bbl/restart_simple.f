c**********************************************************************
      subroutine restart(nhist,to)                                         
c**********************************************************************
c                                                                       
      include 'dim.h'
c
      include 'comwave.h'
      include 'comleng.h'
      include 'comtime.h'
      include 'comflow.h'
      include 'comparamphys.h'
      include 'comcoeff.h'
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
!      include 'comstratif.h' ! TAK ADDED 
   
C
      read(nhist) to
      read(nhist) pr,ra
      read(nhist) li,lj,lk,lkm,lippj,lijk,lsubd,lzloc,lsdtype
c
      write(*,*) 'READING IN RESTART FILE'
c
c
      if( (li.ne.nx) .or. (lj.ne.ny) .or. (lk.ne.nz) ) then
         print *,'warning !!!!'
         print *,'dimensions are inconsistent'
         print *,'li,lj,lk: ',li,lj,lk
         print *,'nx,ny,nz: ',nx,ny,nz
      endif
c
      read(nhist) alpha,beta                                            
      read(nhist) (xw(i),xsq(i),i=1,li)                             
      read(nhist) (yw(i),ysq(i),i=1,lj)                             
      read(nhist) (x(i),i=1,nx)
      read(nhist) (y(j),j=1,ny)
      read(nhist) (z(k),k=1,nz)
      read(nhist) (zf(k),k=1,nz2)
      read(nhist) (z0(k),zh(k),isubdtype(k),k=1,nsubd)
      read(nhist) (subdlen(k),k=1,nsdtype) 
      read(nhist) ((alpol(j,k),j=0,nzm),k=0,nzm)
      read(nhist) ((dalpol(j,k),j=0,nzm),k=0,nzm)
      read(nhist) (zpts(i),i=0,lkm)                                     
      read(nhist) (wg(i),i=0,lkm)
      read(nhist) ((d(i,j),i=0,lkm),j=0,lkm)
      read(nhist) ((d2(i,j),i=0,lkm),j=0,lkm)
      read(nhist) ((d3(i,j),i=0,lkm),j=0,lkm)
      read(nhist) (atpz(k),btpz(k),ctpz(k),k=1,nz)
      read(nhist) ((tubc(i,j),tbbc(i,j),i=1,nxpp),j=1,ny) 
      read(nhist) (((uw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2),
     &             (((vw(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2),
     &             (((ww(i,j,k),i=1,2*nxpp),j=1,ny),k=1,2) 

c Soliton velocity/density and gradient fields (includes pressure gradient)
      read(nhist) umax ! TAK ADDED
      read(nhist) solks
      read(nhist) ((umwv(i,k),wmwv(i,k),i=1,nxpp),k=1,nz)
      read(nhist) ((dumdx(i,k),dumdz(i,k),i=1,nxpp),k=1,nz)
      read(nhist) ((dwmdx(i,k),dwmdz(i,k),i=1,nxpp),k=1,nz)
      read(nhist) ((rhomwv(i,k),i=1,nxpp),k=1,nz)
      read(nhist) ((drhomdx(i,k),drhomdz(i,k),i=1,nxpp),k=1,nz)
!     read(nhist) (rhobar(k),k=1,nz)      ! TAK ADDED
!     read(nhist) (drhobardz(k),k=1,nz)   ! TAK ADDED
!     read(nhist) (d2rhobardz2(k),k=1,nz) ! TAK ADDED
      
c
c velocity and temperature at (n) time level in spectral-physical space
c
      read(nhist) ( ((u(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((v(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((w(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((temp(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)

c
c velocity and temperature (n-1) time level in physical space
c (for BDF purposes)
c
      read(nhist) ( ((un(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((vn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((wn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((tn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)

c
c velocity and temperature (n-2) time level in physical space
c (for BDF purposes)
c
      read(nhist) ( ((unm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((vnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((wnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((tnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
 

c
c convective terms at (n-1) time level in physical space
c
      read(nhist) ( ((sun(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((svn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((swn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((stn(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)

c
c convective terms at (n-1) time level in physical space
c
      read(nhist) ( ((sunm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((svnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((swnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
      read(nhist) ( ((stnm(i,j,k),i=1,nxpp),j=1,ny),k=1,nz)
c
 
c     do k =1,nz
c       i = nxh
c       j = nyh+1  
c        write(*,'(5(E12.6,2x))') 
c    >  z(k),u(i,j,k),v(i,j,k),w(i,j,k),temp(i,j,k)
c     enddo

      return                                                            
      end                                                               
