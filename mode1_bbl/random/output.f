!**********************************************************************
      subroutine output(t,fout,myid,id1d,id2d,comm,
     >                  uf,vf,wf,tempf,
     >                  un,vn,wn,tn,
     >                  unm,vnm,wnm,tnm,
     >                  sun,svn,swn,stn,
     >                  sunm,svnm,swnm,stnm,istep)
!**********************************************************************
!-Simplified version, where each processor dumps out its own file 
!-containing specific block of data.

#if USE_NETCDF
      use netcdf
#endif

#if USE_MPI
      use mpi
#endif

      include 'dim.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comparamphys.h'
      include 'combc.h'
      include 'comlegendre.h'
      include 'comsubd.h'
      include 'comgridptxyz.h'
      include 'comspecfilter.h'
      include 'comsubdpar.h'
      include 'comverbosity.h'
      integer :: verbosity
!-3-D Arrays (GLOBAL)
!-Fourier Space only !
      real, dimension(nxpl,ny,nzpl) :: uf,vf,wf,tempf,
     >                                 un,vn,wn,tn,
     >                                 unm,vnm,wnm,tnm,
     >                                 sun,svn,swn,stn,
     >                                 sunm,svnm,swnm,stnm

#if USE_NETCDF
!- Define variables needed by NETCDF output routines -- Kris R.
      
      integer :: ncid, FLAG, chunksize
      integer :: starts(3), counts(3), chunksizes(3) !This should get changed if we want to make code 2D

!- File name
      character (len=8) :: plot_index
      character (len=20) :: FILE_NAME

!- Dimension names
      character (len=*), parameter :: X_NAME = "x"
      character (len=*), parameter :: Y_NAME = "y"
      character (len=*), parameter :: Z_NAME = "z"
      integer :: x_dimid, y_dimid, z_dimid
      integer :: dimids(3) !This should get changed if we want to make code 2D

!- Data field names
      character (len=*), parameter :: ITER = "Cycle"
      character (len=*), parameter :: TIME = "Time"
      character (len=*), parameter :: XLENGTH = "xlen"
      character (len=*), parameter :: YLENGTH = "ylen"
      character (len=*), parameter :: ZLENGTH = "zlen"
      character (len=*), parameter :: NUMSUBD = "nsubd"
      character (len=*), parameter :: NZLOCATIONS = "nzloc"
      
      character (len=*), parameter :: UF_NAME = "uf"
      character (len=*), parameter :: VF_NAME = "vf"
      character (len=*), parameter :: WF_NAME = "wf"
      character (len=*), parameter :: TF_NAME = "tf"
      character (len=*), parameter :: UN_NAME = "un"
      character (len=*), parameter :: VN_NAME = "vn"
      character (len=*), parameter :: WN_NAME = "wn"
      character (len=*), parameter :: TN_NAME = "tn"
      character (len=*), parameter :: UNM_NAME = "unm"
      character (len=*), parameter :: VNM_NAME = "vnm"
      character (len=*), parameter :: WNM_NAME = "wnm"
      character (len=*), parameter :: TNM_NAME = "tnm"
      character (len=*), parameter :: SUN_NAME = "sun"
      character (len=*), parameter :: SVN_NAME = "svn"
      character (len=*), parameter :: SWN_NAME = "swn"
      character (len=*), parameter :: STN_NAME = "stn"
      character (len=*), parameter :: SUNM_NAME = "sunm"
      character (len=*), parameter :: SVNM_NAME = "svnm"
      character (len=*), parameter :: SWNM_NAME = "swnm"
      character (len=*), parameter :: STNM_NAME = "tvnm"
      integer :: uf_varid, vf_varid, wf_varid, tf_varid
      integer :: un_varid, vn_varid, wn_varid, tn_varid
      integer :: unm_varid, vnm_varid, wnm_varid, tnm_varid
      integer :: sun_varid, svn_varid, swn_varid, stn_varid
      integer :: sunm_varid, svnm_varid, swnm_varid, stnm_varid

      character (len=*), parameter :: PROCID = "mpi_rank"
      integer :: x_varid, y_varid, z_varid, rank_varid

      real, dimension(nxpl,ny,nzpl) :: krisOut
#else
      character (len=39) :: fout
      character (len=8) :: fmt1
#endif

!-------------
!MPI Variables
!-------------
      integer :: myid, ierr, comm
      integer, dimension(nproch*nprocv,2) :: id2d
      integer, dimension(nproch,nprocv) :: id1d

      if((myid.eq.0).and.(verbosity.gt.1)) then
         print *, "----------------------------------------"
         print *, "  In subroutine OUTPUT"
      endif !verbosity > 1

!- If we are using NETCDF execute this code
#if USE_NETCDF
      write(plot_index, "(I6.6)") istep
      FILE_NAME = "chk."//trim(plot_index)//".nc"
!- Create the netcdf file 

!- Set the chunk size. This is very important for parallel I/O performance.
!- This should get revisited and tested more thoroughly at some point
      chunksizes = (/ nxpl, ny, nzpl/) 
#if USE_MPI
      FLAG = IOR(nf90_clobber, NF90_NETCDF4)
      FLAG = IOR(NF90_MPIIO, FLAG)
      call check(nf90_create(FILE_NAME, FLAG, ncid,
     &                  comm=MPI_COMM_WORLD,info=MPI_INFO_NULL))
#else
      print *, "NON-MPI NETCDF NOT IMPLEMENTED YET!!"
      stop "Stopped!"
#endif

!- Define the dimensions
      call check(nf90_def_dim(ncid, X_NAME, nx, x_dimid))   
      call check(nf90_def_dim(ncid, Y_NAME, ny, y_dimid))   
      call check(nf90_def_dim(ncid, Z_NAME, nzpl, z_dimid))   

!- Define the coordinate variables
      call check(nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, x_varid))
      call check(nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, y_varid))
      call check(nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, z_varid))


!- Define the data fields
      dimids = (/x_dimid,y_dimid,z_dimid/)
#if 0
      call check(nf90_def_var(ncid, PROCID, NF90_INT, dimids, rank_varid))

      call check( nf90_def_var(ncid, UF_NAME, NF90_REAL, dimids, uf_varid))
      call check( nf90_def_var(ncid, VF_NAME, NF90_REAL, dimids, vf_varid))
      call check( nf90_def_var(ncid, WF_NAME, NF90_REAL, dimids, wf_varid))
      call check( nf90_def_var(ncid, TF_NAME, NF90_REAL, dimids, tf_varid))

      call check( nf90_def_var(ncid, UN_NAME, NF90_REAL, dimids, un_varid))
      call check( nf90_def_var(ncid, VN_NAME, NF90_REAL, dimids, vn_varid))
      call check( nf90_def_var(ncid, WN_NAME, NF90_REAL, dimids, wn_varid))
      call check( nf90_def_var(ncid, TN_NAME, NF90_REAL, dimids, tn_varid))
      
      call check( nf90_def_var(ncid, UNM_NAME, NF90_REAL, dimids, unm_varid))
      call check( nf90_def_var(ncid, VNM_NAME, NF90_REAL, dimids, vnm_varid))
      call check( nf90_def_var(ncid, WNM_NAME, NF90_REAL, dimids, wnm_varid))
      call check( nf90_def_var(ncid, TNM_NAME, NF90_REAL, dimids, tnm_varid))
      
      call check( nf90_def_var(ncid, SUN_NAME, NF90_REAL, dimids, sun_varid))
      call check( nf90_def_var(ncid, SVN_NAME, NF90_REAL, dimids, svn_varid))
      call check( nf90_def_var(ncid, SWN_NAME, NF90_REAL, dimids, swn_varid))
      call check( nf90_def_var(ncid, STN_NAME, NF90_REAL, dimids, stn_varid))
      
      call check( nf90_def_var(ncid, SUNM_NAME, NF90_REAL, dimids, sunm_varid))
      call check( nf90_def_var(ncid, SVNM_NAME, NF90_REAL, dimids, svnm_varid))
      call check( nf90_def_var(ncid, SWNM_NAME, NF90_REAL, dimids, swnm_varid))
      call check( nf90_def_var(ncid, STNM_NAME, NF90_REAL, dimids, stnm_varid))
#endif     
 
      call check(nf90_def_var(ncid, PROCID, NF90_INT, dimids, 
     &            rank_varid, chunksizes=chunksizes))

      call check( nf90_def_var(ncid, UF_NAME, NF90_REAL, dimids, 
     &            uf_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, VF_NAME, NF90_REAL, dimids, 
     &            vf_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, WF_NAME, NF90_REAL, dimids, 
     &            wf_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, TF_NAME, NF90_REAL, dimids, 
     &            tf_varid, chunksizes=chunksizes))

      call check( nf90_def_var(ncid, UN_NAME, NF90_REAL, dimids, 
     &            un_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, VN_NAME, NF90_REAL, dimids, 
     &            vn_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, WN_NAME, NF90_REAL, dimids, 
     &            wn_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, TN_NAME, NF90_REAL, dimids, 
     &            tn_varid, chunksizes=chunksizes))
      
      call check( nf90_def_var(ncid, UNM_NAME, NF90_REAL, dimids, 
     &            unm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, VNM_NAME, NF90_REAL, dimids, 
     &            vnm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, WNM_NAME, NF90_REAL, dimids, 
     &            wnm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, TNM_NAME, NF90_REAL, dimids, 
     &            tnm_varid, chunksizes=chunksizes))
      
      call check( nf90_def_var(ncid, SUN_NAME, NF90_REAL, dimids, 
     &            sun_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, SVN_NAME, NF90_REAL, dimids, 
     &            svn_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, SWN_NAME, NF90_REAL, dimids, 
     &            swn_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, STN_NAME, NF90_REAL, dimids, 
     &            stn_varid, chunksizes=chunksizes))
      
      call check( nf90_def_var(ncid, SUNM_NAME, NF90_REAL, dimids, 
     &            sunm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, SVNM_NAME, NF90_REAL, dimids, 
     &            svnm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, SWNM_NAME, NF90_REAL, dimids, 
     &            swnm_varid, chunksizes=chunksizes))
      call check( nf90_def_var(ncid, STNM_NAME, NF90_REAL, dimids, 
     &            stnm_varid, chunksizes=chunksizes))

!- Write smaller metadata first
      call check( nf90_put_att(ncid, NF90_GLOBAL, ITER,istep))
      call check( nf90_put_att(ncid, NF90_GLOBAL, TIME,t))
      call check( nf90_put_att(ncid, NF90_GLOBAL, XLENGTH,xlen))
      call check( nf90_put_att(ncid, NF90_GLOBAL, YLENGTH,ylen))
      call check( nf90_put_att(ncid, NF90_GLOBAL, ZLENGTH,zlen))
      call check( nf90_put_att(ncid, NF90_GLOBAL, NUMSUBD,nsubd))
      call check( nf90_put_att(ncid, NF90_GLOBAL, NZLOCATIONS,nzloc))
      
      call check(nf90_put_var(ncid, x_varid, x))
      call check(nf90_put_var(ncid, y_varid, yglob))
      call check(nf90_put_var(ncid, z_varid, z))
      

!- Since each processor contains one slab of data in the y-direction we
!- need to tell NETCDF to stride the write operations so that everything
!- ends up in the correct place. -- Kris R.

      starts = (/ myid*nxpl+1, 1, 1/)
      counts = (/nxpl, ny, nzpl/)
      
      krisOut(1:nxpl,1:ny,1:nzpl) = myid
      call check(nf90_put_var(ncid, rank_varid, krisOut, starts, counts))

!-  Write velocity and temperature at t^n in spectral space     
      call check(nf90_put_var(ncid, uf_varid, uf, starts, counts))
      call check(nf90_put_var(ncid, vf_varid, vf, starts, counts))
      call check(nf90_put_var(ncid, wf_varid, wf, starts, counts))
      call check(nf90_put_var(ncid, tf_varid, tempf, starts, counts))
      

!-  Write velocity and temperature at t^(n-1) in physical space     
      call check(nf90_put_var(ncid, un_varid, un, starts, counts))
      call check(nf90_put_var(ncid, vn_varid, vn, starts, counts))
      call check(nf90_put_var(ncid, wn_varid, wn, starts, counts))
      call check(nf90_put_var(ncid, tn_varid, tn, starts, counts))

!-  Write velocity and temperature at t^(n-2) in physical space     
      call check(nf90_put_var(ncid, unm_varid, unm, starts, counts))
      call check(nf90_put_var(ncid, vnm_varid, vnm, starts, counts))
      call check(nf90_put_var(ncid, wnm_varid, wnm, starts, counts))
      call check(nf90_put_var(ncid, tnm_varid, tnm, starts, counts))

!-  Write advection terms at t^(n-1) in physical space     
      call check(nf90_put_var(ncid, sun_varid, sun, starts, counts))
      call check(nf90_put_var(ncid, svn_varid, svn, starts, counts))
      call check(nf90_put_var(ncid, swn_varid, swn, starts, counts))
      call check(nf90_put_var(ncid, stn_varid, stn, starts, counts))

!-  Write advection at t^(n-2)
      call check(nf90_put_var(ncid, sunm_varid, sunm, starts, counts))
      call check(nf90_put_var(ncid, svnm_varid, svnm, starts, counts))
      call check(nf90_put_var(ncid, swnm_varid, swnm, starts, counts))
      call check(nf90_put_var(ncid, stnm_varid, stnm, starts, counts))

!- Finally close the file 
      call check(nf90_close(ncid))

!- If we are not using NETCDF then do things the old way.
#else
!-Create output file
!-CAREFUL: Naming scheme can handle up to 999 processors !
      if ((myid >= 0) .and. (myid < 10) )    fmt1 = '(a3,i1)'
      if ((myid >= 10) .and. (myid < 100))   fmt1 = '(a2,i2)'
      if ((myid >= 100) .and. (myid < 1000)) fmt1 = '(a1,i3)'
  
      if ((myid >= 0) .and. (myid < 10) )  
     >write(unit=fout(36:39),fmt=fmt1) '_00',myid
      if ((myid >= 10) .and. (myid < 100))
     >write(unit=fout(36:39),fmt=fmt1) '_0',myid
      if ((myid >= 100) .and. (myid < 1000))
     >write(unit=fout(36:39),fmt=fmt1) '_',myid

!-Open output file (each processor does that separately, on his own channel)
      iunit = 2000+myid
      open(iunit,file=fout,form='unformatted',status='unknown')
!     write(*,*) 'proc & file',myid,fout

!-Write to output file
      write(iunit) t,xlen,ylen,zlen
      write(iunit) nx,ny,nz,nzloc,nsubd,nsdtype,nproch,nprocv
     

! velocity and temperature at (n) time level in spectral-physical space
!
      write(iunit) ( ((uf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((vf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((wf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((tempf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
! velocity and temperature (n-1) time level in physical space
! (for BDF purposes)
!
      write(iunit) ( ((un(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((vn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((wn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((tn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
!
! velocity and temperature (n-2) time level in physical space
! (for BDF purposes)
!
      write(iunit) ( ((unm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((vnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((wnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((tnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
! convective terms at (n-1) time level in physical space
!
      write(iunit) ( ((sun(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((svn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((swn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((stn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!
!
! convective terms at (n-1) time level in physical space
!
      write(iunit) ( ((sunm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((svnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((swnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      write(iunit) ( ((stnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!-Close output file
      close(iunit)

!      call MPI_BARRIER(comm,ierr)

!- End if not using NETCDF
#endif

      end subroutine output 

!**********************************************************************
      subroutine plotf(t,fout,myid,numprocs,comm,twoslice,subslice,
     >                 id1d,id2d,
     >                 u,v,w,temp,uf,vf,wf,tempf,plot_prefix,istep)
!**********************************************************************
!
!-Subroutine that dumps out 3-D fields of (u,v,w,rho) in
!-SINGLE PRECISION for future postprocessing. Uses MPI-2
!-in exactly the same fashion as subroutine *output* does.
!-Created by PD on 8/8/05.
!-NOTE: As in serial version of code data is outputted in
!-Physical Space. Though this may entail some additional FFTs
!-it makes the data more "shareable" with the community.

#if USE_NETCDF
      use netcdf
!      use plot
#endif

#if USE_MPI
      use mpi
#endif

      include 'dim.h'
      include 'comwave.h'
      include 'comleng.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comsubdpar.h'
!      include 'comtime.h'
      include 'comparamphys.h'
      include 'comverbosity.h'
      integer :: verbosity

!-3-D Arrays (GLOBAL)
!-a) Phys. Space
      real, dimension(nxpp,nypl,nzpl) :: u,v,w,temp
!-b) Fourier Space
      real, dimension(nxpl,ny,nzpl) :: uf,vf,wf,tempf

#if USE_NETCDF
!- Define variables needed by NETCDF output routines -- Kris R.
      
      integer :: ncid, FLAG
      integer :: starts(3), counts(3) !This should get changed if we want to make code 2D

!- File name
      character (len=8) :: plot_index
      character (len=40) :: FILE_NAME
      character (len=25) :: plot_prefix

!- Dimension names
      character (len=*), parameter :: X_NAME = "x"
      character (len=*), parameter :: Y_NAME = "y"
      character (len=*), parameter :: Z_NAME = "z"
      integer :: x_dimid, y_dimid, z_dimid
      integer :: dimids(3) !This should get changed if we want to make code 2D

!- Coordinate fields
      integer :: x_varid, y_varid, z_varid

!- Data field names
      character (len=*), parameter :: ITER = "Cycle"
      character (len=*), parameter :: TIME = "Time"
      character (len=*), parameter :: NT = "Nt"
      character (len=*), parameter :: D_TIME = "dt"
      character (len=*), parameter :: U_NAME = "x-velocity"
      character (len=*), parameter :: V_NAME = "y-velocity"
      character (len=*), parameter :: W_NAME = "z-velocity"
      character (len=*), parameter :: T_NAME = "temperature"
      integer :: u_varid, v_varid, w_varid, t_varid
      integer :: rank_varid, subd_varid

      character (len=*), parameter :: PROCID = "mpi_rank"
      character (len=*), parameter :: SUBDOMAIN = "subdomain"
      integer :: subdindx, subdid

!- Include units for dimensions and data fields
      character (len=*), parameter :: UNITS  = "units"
      character (len=*), parameter :: X_UNITS="m"
      character (len=*), parameter :: Y_UNITS="m"
      character (len=*), parameter :: Z_UNITS="m"
      character (len=*), parameter :: U_UNITS="ms^-1"
      character (len=*), parameter :: V_UNITS="ms^-1"
      character (len=*), parameter :: W_UNITS="ms^-1"
      character (len=*), parameter :: T_UNITS="K"

#else
      character (len=60) :: fout
!      character (len=42) :: fout
      character (len=8) :: fmt1
#endif

!-------------
!MPI Variables
!-------------
      integer :: myid, numprocs, ierr, comm
      integer :: twoslice, subslice
      integer, dimension(nproch*nprocv,2) :: id2d
      integer, dimension(nproch,nprocv) :: id1d

!-Temporary output array
!-NOTE: Kind = 4 is compiler dependent !!
!      real(kind=4), dimension(:,:,:), allocatable :: tempout
      real, dimension(:,:,:), allocatable :: tempout
      real :: krisOut(nx,nypl,nzpl)

      integer AllocateStatus

      if((myid.eq.0).and.(verbosity.gt.1)) then
         print *, "----------------------------------------"
         print *, "  In subroutine PLOTF"
      endif !verbosity > 1

!-Allocate local arrays
      allocate ( tempout(nxpp,nypl,nzpl),
     >           stat = AllocateStatus)
                                                                                
      if (AllocateStatus /= 0) then
         stop "**Not Enough Memory - OUTPUT**"
      end if

!-Before doing anything else convert flow field back to
!-Physical Space
      call horfft(u,uf,1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(v,vf,1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(w,wf,1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(temp,tempf,1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)

#if USE_NETCDF
      write(plot_index, "(I6.6)") istep
      FILE_NAME = trim(plot_prefix)//"."//trim(plot_index)//".nc"
      FILE_NAME = trim(FILE_NAME)

!- Create the netcdf file
#if USE_MPI
      FLAG = IOR(nf90_clobber,NF90_NETCDF4)
      FLAG = IOR(NF90_MPIIO, FLAG)
      call check(nf90_create(FILE_NAME, FLAG, ncid,  
     &               comm=comm,info=MPI_INFO_NULL))
     
#else
      print *, "NON-MPI NETCDF NOT IMPLEMENTED YET!!"
      stop "Stopped!"
#endif

!- Define the dimensions
      call check(nf90_def_dim(ncid, X_NAME, nx, x_dimid))   
      call check(nf90_def_dim(ncid, Y_NAME, ny, y_dimid))   
      call check(nf90_def_dim(ncid, Z_NAME, nzpl, z_dimid))   

!- Define the coordinate variables
      call check(nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, x_varid))
      call check(nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, y_varid))
      call check(nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, z_varid))

!- Assign units to each coordinate variable
      call check(nf90_put_att(ncid, x_varid, UNITS, X_UNITS))
      call check(nf90_put_att(ncid, y_varid, UNITS, Y_UNITS))
      call check(nf90_put_att(ncid, z_varid, UNITS, Z_UNITS))

!- Define the data fields
      dimids = (/x_dimid,y_dimid,z_dimid/)
      call check( nf90_def_var(ncid, U_NAME, NF90_REAL, dimids, u_varid))
      call check( nf90_def_var(ncid, V_NAME, NF90_REAL, dimids, v_varid))
      call check( nf90_def_var(ncid, W_NAME, NF90_REAL, dimids, w_varid))
      call check( nf90_def_var(ncid, T_NAME, NF90_REAL, dimids, t_varid))

      call check(nf90_def_var(ncid, PROCID, NF90_INT, dimids, rank_varid))
      call check(nf90_def_var(ncid, SUBDOMAIN, NF90_INT, dimids, subd_varid))

!- Assign units to each of the data fields

      call check(nf90_put_att(ncid, u_varid, UNITS, U_UNITS))
      call check(nf90_put_att(ncid, v_varid, UNITS, V_UNITS))
      call check(nf90_put_att(ncid, w_varid, UNITS, w_UNITS))
      call check(nf90_put_att(ncid, t_varid, UNITS, T_UNITS))
      
!- Write the time and cycle number as global attributes
      call check( nf90_put_att(ncid, NF90_GLOBAL, ITER,istep))
      call check( nf90_put_att(ncid, NF90_GLOBAL, TIME,t))
      call check( nf90_put_att(ncid, NF90_GLOBAL, NT, 
     & brunt*(t-tstart)+brunt0))
      call check( nf90_put_att(ncid, NF90_GLOBAL, D_TIME,dt))

!- Since each processor contains one slab of data in the y-direction we
!- need to tell NETCDF to stride the write operations so that everything
!- ends up in the correct place. -- Kris R.

      starts = (/ 1, myid*nypl+1 , 1/)
      counts = (/nx, nypl, nzpl/)

!- Now write the coordinate variables
      call check(nf90_put_var(ncid, x_varid, x))
      call check(nf90_put_var(ncid, y_varid, yglob))
      call check(nf90_put_var(ncid, z_varid, z))

!- Write out the processor IDs
      
      krisOut(1:nx, 1:nypl, 1:nzpl) = myid
      call check(nf90_put_var(ncid, rank_varid, krisOut, starts, counts))

!- Write out the subdomain IDs
         do k=1,nzpl
               
               subdid = 0  
               do subdindx = 1, nsubd                             
                  if (z(k).ge.z0(subdindx)) then
                     subdid = subdid + 1
                  endif
               enddo
      
               krisOut(1:nx,1:nypl,k) = subdid
                                            
         enddo

      call check(nf90_put_var(ncid, subd_varid, krisOut, starts, counts))

!- Write the data fields
!- Kris R. -- For now just do it the quick and dirty way pete has in the code
!- already, Afterwards do something that avoids extra FFTs at the end.


!-Dump u-velocity
!      call horfft(u,uf,1,
!     >            myid,numprocs,comm,
!     >            twoslice,subslice,id1d,id2d)
      krisOut(1:nx,1:nypl,1:nzpl)=real(u(2:(nx+1),1:nypl,1:nzpl))
      call check(nf90_put_var(ncid, u_varid, krisOut, starts, counts))

!-Dump v-velocity
!      call horfft(v,vf,1,
!     >            myid,numprocs,comm,
!     >            twoslice,subslice,id1d,id2d)
      krisOut(1:nx,1:nypl,1:nzpl)=real(v(2:(nx+1),1:nypl,1:nzpl))
      call check(nf90_put_var(ncid, v_varid, krisOut, starts, counts))

!-Dump w-velocity
!      call horfft(w,wf,1,
!     >            myid,numprocs,comm,
!     >            twoslice,subslice,id1d,id2d)
      krisOut(1:nx,1:nypl,1:nzpl)=real(w(2:(nx+1),1:nypl,1:nzpl))
      call check(nf90_put_var(ncid, w_varid, krisOut, starts, counts))

!-Dump temperature
!      call horfft(temp,tempf,1,
!     >            myid,numprocs,comm,
!     >            twoslice,subslice,id1d,id2d)
      krisOut(1:nx,1:nypl,1:nzpl)=real(temp(2:(nx+1),1:nypl,1:nzpl))
      call check(nf90_put_var(ncid, t_varid, krisOut, starts, counts))

!- Finally close the file 
      call check(nf90_close(ncid))

#else 

!-Update animation file list
      if (myid == 0) write(69,'(a10,a12)') 'postp_dir/',fout(10:21)
!-Create output file
!-CAREFUL: Naming scheme can handle up to 999 processors !
      if ((myid >= 0) .and. (myid < 10) )    fmt1 = '(a3,i1)'
      if ((myid >= 10) .and. (myid < 100))   fmt1 = '(a2,i2)'
      if ((myid >= 100) .and. (myid < 1000)) fmt1 = '(a1,i3)'

      if ((myid >= 0) .and. (myid < 10) )
     >write(unit=fout(39:42),fmt=fmt1) '_00',myid
      if ((myid >= 10) .and. (myid < 100))
     >write(unit=fout(39:42),fmt=fmt1) '_0',myid
      if ((myid >= 100) .and. (myid < 1000))
     >write(unit=fout(39:42),fmt=fmt1) '_',myid

!-Open output file (each processor does that separately, on his own channel)
      iunit = 2000+myid
      open(iunit,file=fout,form='unformatted',status='unknown')

!-First dump out information on grid resolution/number of processors
!-and dump out some physical constants etc.
!--------------------------------------
       write(iunit) t,xlen,ylen,zlen
       write(iunit) nx,ny,nz,nzloc,nsubd,nsdtype,nproch,nprocv

!-NOW DUMP OUT 3-D FIELDS

!-The temporary copies are not done through some subroutine
!-to avoid screw-ups with memory addressing (PD: Improve this asap !)
!-Dump u-velocity
         do k=1,nzpl
           do j=1,nypl
             do i=1,nxpp
                                                                                
               tempout(i,j,k) = real(u(i,j,k))
                                                                                
             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,nypl),k=1,nzpl)

!-Dump v-velocity
         do k=1,nzpl
           do j=1,nypl
             do i=1,nxpp
                                                                                
               tempout(i,j,k) = real(v(i,j,k))
                                                                                
             enddo
           enddo
         enddo

         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,nypl),k=1,nzpl)      

!-Dump w-velocity
         do k=1,nzpl
           do j=1,nypl
             do i=1,nxpp
                                                                                
               tempout(i,j,k) = real(w(i,j,k))
                                                                                
             enddo
           enddo
         enddo
                                                                                
        write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,nypl),k=1,nzpl)

!-Dump Active Scalar
         do k=1,nzpl
           do j=1,nypl
             do i=1,nxpp
                                                                                
               tempout(i,j,k) = real(temp(i,j,k))
                                                                                
             enddo
           enddo
         enddo
                                                                                
         write(iunit) (((tempout(i,j,k),i=1,nxpp),j=1,nypl),k=1,nzpl)


!-Close output file
         close(iunit)

#endif

!         call MPI_BARRIER(comm,ierr)


!-Convert back to Fourier space
      call horfft(u,uf,-1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(v,vf,-1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(w,wf,-1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
      call horfft(temp,tempf,-1,
     >            myid,numprocs,comm,
     >            twoslice,subslice,id1d,id2d)
 
!-De-Allocate local array
      deallocate (tempout,
     >            stat=AllocateStatus)
      if (AllocateStatus /= 0) then
         stop "**Error Deallocating - OUTPUT**"
      end if
  
      end subroutine plotf

#if USE_NETCDF

!-    A subroutine to check error flags returned from NETCDF routines
      subroutine check(status)
         use netcdf

         integer, intent (in) :: status

         if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped!"
         endif

      end subroutine check

#endif

