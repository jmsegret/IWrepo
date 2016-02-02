c**********************************************************************
      subroutine restart(nhist,to)                                         
c**********************************************************************
c                                                                       
      include 'dim.h'
c
      include 'comtime.h'
      include 'comwave.h'
      include 'comleng.h'
c      
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

             parameter(Nxpold=nxpp)
!              parameter(Nxpold=514)

c     write(*,*)'---------------------------'
c     write(*,*)'Regridding with nx=',Nxpold-2

c nxpp=nx+2=128+2=130
C if Regridding ie performed Nxpold is the old nx+2
C*** A.M.A. Oct-13-08 
        real, dimension(:,:,:), allocatable :: temp3drg  
        integer AllocateStatus

       common /Dump/ Tdump,Tx        
       common /anim/ianim,ianimc 
C
        allocate(temp3drg(Nxpold,ny,nz),stat=AllocateStatus)
       if (AllocateStatus /= 0) then
         stop "**Not Enough Memory -Horizontal Regridding**"
       endif
      write(*,*)'----------------------------------------------'
      write(*,*)'Regridding with nx=',nx,'From nx=',Nxpold-2
      write(*,*)'----------------------------------------------'



      read(nhist) to
      read(nhist) pr,ra
      read(nhist) li,lj,lk,lkm,lippj,lijk,lsubd,lzloc,lsdtype
c
      write(*,*) 'READING IN RESTART FILE'
c
c
      if( (li.ne.Nxpold-2) .or. (lj.ne.ny) .or. (lk.ne.nz) ) then
         print *,'warning !!!!'
         print *,'dimensions are inconsistent'
         print *,'li,lj,lk: ',li,lj,lk
         print *,'nx,ny,nz: ',Nxpold-2,ny,nz
      endif

      read(nhist) alpha,beta                                            
      read(nhist) (xw(i),xsq(i),i=1,li)                             
      read(nhist) (yw(i),ysq(i),i=1,lj)       
      read(nhist) (x(i),i=1,Nxpold-2)
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
C************************************************************************************
       read(nhist) ((temp3drg(i,j,1),i=1,Nxpold),j=1,ny)
       call regcopy2d(temp3drg,Nxpold,nxpp,ny,nz,tubc)     
       read(nhist) ((temp3drg(i,j,1),i=1,Nxpold),j=1,ny)      
       call regcopy2d(temp3drg,Nxpold,nxpp,ny,nz,tbbc)
       read(nhist) (((temp3drg(i,j,k),i=1,2*Nxpold),j=1,ny),k=1,2)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,2,uw)
       read(nhist)(((temp3drg(i,j,k),i=1,2*Nxpold),j=1,ny),k=1,2)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,2,vw)
       read(nhist)(((temp3drg(i,j,k),i=1,2*Nxpold),j=1,ny),k=1,2)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,2,ww)
c velocity and temperature at (n) time level in spectral-physical space
       read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,u)
       read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,v)
       read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,w)
       read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
       call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,temp)

c velocity and temperature (n-1) time level in physical space
c (for BDF purposes)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,un)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,vn)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,wn)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,tn)
c velocity and temperature (n-2) time level in physical space
c (for BDF purposes)
c
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,unm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,vnm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,wnm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,tnm)

c convective terms at (n-1) time level in physical space
c
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,sun)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,svn)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,swn)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,stn)
c convective terms at (n-1) time level in physical space
c
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,sunm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,svnm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,swnm)
      read(nhist) ( ((temp3drg(i,j,k),i=1,Nxpold),j=1,ny),k=1,nz)
      call regcopy3d(temp3drg,Nxpold,nxpp,ny,nz,stnm)

      read(nhist) ianimc,Tdump,Tx
C**************************************
       deallocate(temp3drg,stat=AllocateStatus)

      if (AllocateStatus /= 0) then
         stop "**Error Deallocating - Horizontal regridding **"
      end if
C*******************************************
      return                                                            
      end                                                    
          
c************************************************************
      subroutine  regcopy3d(Uinp,Nxpold,nxpp,ny,nz,Uoutp)
      
         real, dimension(1:Nxpold,ny,nz) ::Uinp  
         real,dimension(1,nxpp,ny,nz) ::Uoutp
          Uoutp=0.0
         do k=1,nz
            do j=1,ny  
                do i=1,Nxpold         
                      Uoutp(i,j,k)=Uinp(i,j,k)
                enddo
            enddo
         end do
          Uinp=0.0
         return
          end
 
       subroutine  regcopy2d(Uinp,Nxpold,nxpp,ny,nz,Uoutp)

         real, dimension(1:Nxpold,ny,nz) ::Uinp
         real,dimension(1,nxpp,ny) ::Uoutp
          Uoutp=0.0
         
            do j=1,ny
                do i=1,Nxpold
                      Uoutp(i,j)=Uinp(i,j,1)
                enddo
            end do
          Uinp=0.0
         return
          end


