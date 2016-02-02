C***********************************************************
       subroutine bcreset
C***********************************************************
C-Specialized subroutine (PD: 12/18/03) which resets the
C-boundary conditions due to change in soliton amplitude
C-(PD: 2/25/04: For internal waves, this only necessitates
C-change in u-BC)
       include 'dim.h'
c
       include 'comflow.h'
       include 'comwave.h'
       include 'comleng.h'
       include 'combc.h'
       include 'comscratch.h'
       include 'comparamphys.h'
       include 'comgridptxyz.h'
       include 'comsubd.h'
       include 'comsoliton.h'

C-(Re)set Boundary conditions (Remember, flow field is in Fourier Space
C-at this point !)


      do j=1,ny
         do i=1,nx

            oy(i,j,1) =  - umwv(i,1)
            oy(i,j,nz) =  0. 

         enddo
      enddo

C-Transform into Fourier space and assign to specific BC
C-variables.
      call  horfft (oy,-1,nz)

      do j=1,ny
         do i=1,nxpp

            uw(i,j,2) = oy(i,j,nz)
            uw(i,j,1) = oy(i,j,1)

c           if (j.eq.1)  write(*,*) uw(i,j,2),ww(i,j,2)
         enddo
      enddo 
      
      return
      end
