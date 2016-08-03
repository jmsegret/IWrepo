!**********************************************************************
      subroutine restart(t,fin,myid,id1d,id2d,comm,
     >                   uf,vf,wf,tempf,
     >                   un,vn,wn,tn,
     >                   unm,vnm,wnm,tnm,
     >                   sun,svn,swn,stn,
     >                   sunm,svnm,swnm,stnm,restart_file,ibegin)
!**********************************************************************
!

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


      integer :: nxin, nyin, nzin, nzlocin, nsubdin
#if USE_NETCDF

      integer :: ncid, FLAG
      integer :: starts(3), counts(3)

      character (len=13) :: FILE_NAME
      character (len=13) :: restart_file

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
      integer :: x_varid, y_varid, z_varid
#else
!-Even in the 2-D case, the global array is considered 3-D
!-with 3rd dimension equal to 1.
      integer, parameter :: ndims_io=3
      character (len=39) :: fin
      character (len=8) :: fmt1
#endif

!-------------
!MPI Variables
!-------------
#if USE_MPI
      integer :: myid, ierr, comm
      integer, dimension(nproch*nprocv,2) :: id2d
      integer, dimension(nproch,nprocv) :: id1d
#endif USE_MPI
      
      if((myid.eq.0).and.(verbosity.gt.1)) then
         print *, '----------------------------------------'
         print *, '  In subroutine RESTART'
      endif !verbosity > 1

!- If using NETCDF compile this code ...
#if USE_NETCDF
!- Kris R: Try to read as much meta-data from the NETCDF file as possible to
!- make restarting capabilities as general as possible. For example, read the time, lengths, and dimensions from the file if possible. This will make tin, etc redundant.

!- Open the checkpoint file
      FILE_NAME=trim(restart_file)
#if USE_MPI
      FLAG = IOR(NF90_NOWRITE,NF90_MPIIO)
      call check(nf90_open(FILE_NAME, FLAG, ncid, 
     &            comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))
#else
      print *, "NON-MPI NETCDF NOT IMPLEMENTED YET!"
      stop "Stopped!"
#endif

#if DEBUG 
!- Check the existence of global attributes
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, ITER))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, TIME))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, XLENGTH))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, YLENGTH))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, ZLENGTH))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, NUMSUBD))
      call check(nf90_inquire_attribute(ncid, NF90_GLOBAL, NZLOCATIONS))
#endif

!- Read in global attributes like the cycle number, time, etc.
      call check(nf90_get_att(ncid, NF90_GLOBAL, ITER, ibegin))
      call check(nf90_get_att(ncid, NF90_GLOBAL, TIME, tstart))
      call check(nf90_get_att(ncid, NF90_GLOBAL, XLENGTH, xlen))
      call check(nf90_get_att(ncid, NF90_GLOBAL, YLENGTH, ylen))
      call check(nf90_get_att(ncid, NF90_GLOBAL, ZLENGTH, zlen))
      t=tstart

!- Kris R: Check to make sure the number of subdomains is the same as in the
!- inputs file
      call check(nf90_get_att(ncid, NF90_GLOBAL, NUMSUBD, nsubdin))
      if(nsubdin .ne. nsubdin) then
         stop "The number of subdomains in the restart file is different
     &      than in dim.h"
      endif

!- Kris R: Check to make sure the number of grid points in each subdomain is
!- the same as in the inputs file
      call check(nf90_get_att(ncid, NF90_GLOBAL, NZLOCATIONS, nzlocin))
      if(nzlocin .ne. nzloc) then
         stop "The number of grid points in each subdomain in the
     &      restart file is different than in dim.h"
      endif

!- Get the id's of dimensions based on their names.
!- If they don't exist an error will be thrown
      call check( nf90_inq_dimid(ncid, X_NAME, x_dimid))
      call check( nf90_inq_dimid(ncid, Y_NAME, y_dimid))
      call check( nf90_inq_dimid(ncid, Z_NAME, z_dimid))

!- Now get the length of each dimension
!- For now require that the number of grid points be the same as in
!- inputs.dat until we change the regridding routine to use netcdf
      call check( nf90_inquire_dimension(ncid, x_dimid, len=nxin))
      if(nxin .ne. nx) then
         stop "The number of x-grid points in the restart file is
     &      different than in dim .h"
      endif

      call check( nf90_inquire_dimension(ncid, y_dimid, len=nyin))
      if(nyin .ne. ny) then
         stop "The number of y-grid points in the restart file is
     &      different than in dim.h"
      endif
      
      call check( nf90_inquire_dimension(ncid, z_dimid, len=nzin))
      if(nzin .ne. nz) then
         stop "The number of z-grid points in the restart file is
     &      different than in dim.h"
      endif

!- Next get the id's of the variables based on their names.
!- If they don't exist then an error will be thrown
      call check( nf90_inq_varid(ncid, UF_NAME, uf_varid))
      call check( nf90_inq_varid(ncid, VF_NAME, vf_varid))
      call check( nf90_inq_varid(ncid, WF_NAME, wf_varid))
      call check( nf90_inq_varid(ncid, TF_NAME, tf_varid))
      
      call check( nf90_inq_varid(ncid, UN_NAME, un_varid))
      call check( nf90_inq_varid(ncid, VN_NAME, vn_varid))
      call check( nf90_inq_varid(ncid, WN_NAME, wn_varid))
      call check( nf90_inq_varid(ncid, TN_NAME, tn_varid))

      call check( nf90_inq_varid(ncid, UNM_NAME, unm_varid))
      call check( nf90_inq_varid(ncid, VNM_NAME, vnm_varid))
      call check( nf90_inq_varid(ncid, WNM_NAME, wnm_varid))
      call check( nf90_inq_varid(ncid, TNM_NAME, tnm_varid))

      call check( nf90_inq_varid(ncid, SUN_NAME, sun_varid))
      call check( nf90_inq_varid(ncid, SVN_NAME, svn_varid))
      call check( nf90_inq_varid(ncid, SWN_NAME, swn_varid))
      call check( nf90_inq_varid(ncid, STN_NAME, stn_varid))

      call check( nf90_inq_varid(ncid, SUNM_NAME, sunm_varid))
      call check( nf90_inq_varid(ncid, SVNM_NAME, svnm_varid))
      call check( nf90_inq_varid(ncid, SWNM_NAME, swnm_varid))
      call check( nf90_inq_varid(ncid, STNM_NAME, stnm_varid))

!- Finally, read the data fields from the file. For now assume that we
!- will use the same number of processors, however later we can change
!- this which will affect how the reads are stridded.
      starts = (/myid*nxpl+1, 1, 1/)
      counts = (/nxpl, ny, nzpl/)

      call check( nf90_get_var(ncid, uf_varid, uf, starts, counts))
      call check( nf90_get_var(ncid, vf_varid, vf, starts, counts))
      call check( nf90_get_var(ncid, wf_varid, wf, starts, counts))
      call check( nf90_get_var(ncid, tf_varid, tempf, starts, counts))

      call check( nf90_get_var(ncid, un_varid, un, starts, counts))
      call check( nf90_get_var(ncid, vn_varid, vn, starts, counts))
      call check( nf90_get_var(ncid, wn_varid, wn, starts, counts))
      call check( nf90_get_var(ncid, tn_varid, tn, starts, counts))
      
      call check( nf90_get_var(ncid, unm_varid, unm, starts, counts))
      call check( nf90_get_var(ncid, vnm_varid, vnm, starts, counts))
      call check( nf90_get_var(ncid, wnm_varid, wnm, starts, counts))
      call check( nf90_get_var(ncid, tnm_varid, tnm, starts, counts))
      
      call check( nf90_get_var(ncid, sun_varid, sun, starts, counts))
      call check( nf90_get_var(ncid, svn_varid, svn, starts, counts))
      call check( nf90_get_var(ncid, swn_varid, swn, starts, counts))
      call check( nf90_get_var(ncid, stn_varid, stn, starts, counts))
      
      call check( nf90_get_var(ncid, sunm_varid, sunm, starts, counts))
      call check( nf90_get_var(ncid, svnm_varid, svnm, starts, counts))
      call check( nf90_get_var(ncid, swnm_varid, swnm, starts, counts))
      call check( nf90_get_var(ncid, stnm_varid, stnm, starts, counts))

!- Close the checkpoint file
      call check(nf90_close(ncid))
    
#if USE_MPI 
      call MPI_BARRIER(comm,ierr)
#endif

      if(verbosity.gt.0) then
         if(myid .eq. 0) then 
            print *, "Restart parameters ... "
            print *, "Last cycle: ", ibegin
            print *, "Start time: ", t
            print *, "xlen: ", xlen
            print *, "ylen: ", ylen
            print *, "zlen: ", zlen
            print *, "nxin: ", nxin
            print *, "nyin: ", nyin
            print *, "nzin: ", nzin
            print *, "nsubdin: ", nsubdin
            print *, "nzlocin: ", nzlocin 
         endif
      endif !verbosity > 0

!- If not using NETCDF, then compile this code instead
#else
!-Create output file
!-CAREFUL: Naming scheme can handle up to 999 processors !
      if ((myid >= 0) .and. (myid < 10) )    fmt1 = '(a3,i1)'
      if ((myid >= 10) .and. (myid < 100))   fmt1 = '(a2,i2)'
      if ((myid >= 100) .and. (myid < 1000)) fmt1 = '(a1,i3)'

      if ((myid >= 0) .and. (myid < 10) )
     >write(unit=fin(36:39),fmt=fmt1) '_00',myid
      if ((myid >= 10) .and. (myid < 100))
     >write(unit=fin(36:39),fmt=fmt1) '_0',myid
      if ((myid >= 100) .and. (myid < 1000))
     >write(unit=fin(36:39),fmt=fmt1) '_',myid

!-Open output file (each processor does that separately, on his own channel)
      iunit = 2000+myid
!      open(iunit,file=fin,form='unformatted',status='unknown')
      open(iunit,file=fin,form='unformatted',status='old')

!-Write to output file
      read(iunit) tin,xlenin,ylenin,zlenin
      read(iunit) nxin,nyin,nzin,nzlocin,nsubdin,
     >nsdtypein,nprochin,nprocvin

!--------------------------------------------
!-Set time
      t = tin

      if (myid == 0) then
        write(*,*)
        write(*,*) '---------------------------------------------'
        write(*,*)
      endif
          
      if (myid == 0) write(*,*) 'PROC-0: Restart time is t=',t

      if (myid == 0) write(*,*) 'PROC-0: Input data dimensions of
     >computational domain',xlenin,ylenin,zlenin   

!-Stop if dimensions are inconsistent 
      if ((xlenin /= xlen).or.
     >    (ylenin /= ylen).or.
     >    (zlenin /= zlen)) then
          if (myid == 0) then
          write(*,*) 'Physical dimensions of computational domain are
     >                not consistent !'
          write(*,*) 'Program is stopped. TRY AGAIN !'
          endif
 
        call MPI_FINALIZE(ierr)
        stop
      endif

      if (myid == 0) write(*,*) 'nx, ny, nz in restart file: ',
     > nxin, nyin, nzin

!-Additional reasons for stopping

      if ((nxin /= nx).or.
     >    (nyin /= ny).or.
     >    (nzin /= nz)) then
        if (myid == 0) then
          write(*,*) 'Total resolution in x,y or z direction
     >              not consistent !'
         write(*,*) 'Program is stopped. TRY AGAIN !'
        endif

        call MPI_FINALIZE(ierr)                                                 
        stop
      endif

      if (myid == 0) write(*,*) 'nzloc, nsubd, nsdtype
     >in restart file: ',
     > nzlocin, nsubdin, nsdtypein
                                                                                
      if ((nzlocin /= nzloc).or.
     >    (nsubdin /= nsubd).or.
     >    (nsdtypein /= nsdtype)) then
        if (myid == 0) then
          write(*,*) 'Order of approximation in individual subdomains,
     >              number of subdomains and types of subdomains not
     >              consistent !'
          write(*,*) 'Program is stopped. TRY AGAIN !'
        endif
                                                                                
        call MPI_FINALIZE(ierr)
        stop
      endif

      if (myid == 0) write(*,*) 'These restart files were 
     >generated with ',
     >nprochin,' horizontal procs. and ', nprocvin, ' vertical procs.'

      if ((nproch /= nprochin).or.
     >    (nprocv /= nprocvin)) then 
        if (myid == 0) then
          write(*,*) 'Number of processors in horizontal or vertical 
     >              direction not consistent !'
          write(*,*) 'Program is stopped. TRY AGAIN !'
        endif

        call MPI_FINALIZE(ierr)
        stop
      endif

      if (myid == 0) then
        write(*,*)
        write(*,*) '---------------------------------------------'
        write(*,*)
      endif      

!-------------------------------------------

!-velocity and temperature at (n) time level in spectral-physical space
!
      read(iunit) ( ((uf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((vf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((wf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((tempf(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
! velocity and temperature (n-1) time level in physical space
! (for BDF purposes)
!
      read(iunit) ( ((un(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((vn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((wn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((tn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
!
! velocity and temperature (n-2) time level in physical space
! (for BDF purposes)
!

      read(iunit) ( ((unm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((vnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((wnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((tnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!
! convective terms at (n-1) time level in physical space
!
      read(iunit) ( ((sun(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((svn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((swn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((stn(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
!
!
! convective terms at (n-1) time level in physical space
!
      read(iunit) ( ((sunm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((svnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((swnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)
      read(iunit) ( ((stnm(i,j,k),i=1,nxpl),j=1,ny),k=1,nzpl)

!-Close output file
      close(iunit)

#if USE_MPI
      call MPI_BARRIER(comm,ierr)
#endif

!end if not using NETCDF
#endif

      end subroutine restart
