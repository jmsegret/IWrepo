C-Implements SGS tensor & vector calculations

C ***********************************************************
      subroutine sgscalc_pre
C ***********************************************************
C-This routine simply sets up some pre-filtering preliminaries
C-for the calculation of SGS fluxes etc.

      include 'dim.h'

      include 'comflow.h'
      include 'comscratch.h'
      include 'comsgs.h'
      include 'commeanv.h'
      include 'comfiltertmp.h'
      

C-SGS MOMENTUM FLUXES
C-Diagonal elements first
C-Copy to temporary arrays first
C-These velocity or velocity*temperature products need to be
C-unfiltered

C-Calculate mean temperature. This is a little redundant
C-because calculation will be repeated in postprocessor segment
C-of code.
C-Go to Fourier Space
        call  horfft (temp,-1,nz)

        do k=1,nz
          tempmean(k)=temp(1,1,k)
        enddo

        call horfft(temp,1,nz)

        do i=1,nx
          do j=1,ny

C-All SGS tensor elements are zero at the bottom !
            t11(i,j,1) = 0.
            t22(i,j,1) = 0.
            t33(i,j,1) = 0.
            t12(i,j,1) = 0.
            t13(i,j,1) = 0.
            t23(i,j,1) = 0.
            th1(i,j,1) = 0.
            th2(i,j,1) = 0.
            th3(i,j,1) = 0.

            do k=2,nz

              t11(i,j,k) = u(i,j,k)*u(i,j,k)
              t22(i,j,k) = v(i,j,k)*v(i,j,k)
              t33(i,j,k) = w(i,j,k)*w(i,j,k)
              t12(i,j,k) = u(i,j,k)*v(i,j,k)
              t13(i,j,k) = u(i,j,k)*w(i,j,k)
              t23(i,j,k) = v(i,j,k)*w(i,j,k)

              ox(i,j,k) = temp(i,j,k) - tempmean(k)
              th1(i,j,k) = u(i,j,k)*ox(i,j,k)
              th2(i,j,k) = v(i,j,k)*ox(i,j,k)
              th3(i,j,k) = w(i,j,k)*ox(i,j,k)

            enddo
          enddo
        enddo

C-Proceed with filtering entire vel/temp field 
C-as well as these products and then
C-continue with SGS calculation
      return
      end



C ***********************************************************
      subroutine sgscalc_post
C ***********************************************************

C-This subroutine does the core of SGS flux calculation

      include 'dim.h'

      include 'comflow.h'
      include 'comscratch.h'
      include 'comsgs.h'
      include 'commeanv.h'
      include 'comfiltertmp.h'      

C-First calculate SGS momentum tensor

C-1) t11 element - (u =/ 0 at top)

      iflag=1
      call modecut(iflag,t11,tmps1)

C-2) t22 element - (v =/0 at top)

      iflag=1
      call modecut(iflag,t22,tmps2)

C-3) t33 element - (w =0 at top)

      iflag=0
      call modecut(iflag,t33,tmps3)

C-4) t12 element - (u =/ 0 at top)

      iflag=1
      call modecut(iflag,t12,tmps4)

C-5) t13 element - (w = 0 at top)

      iflag=0
      call modecut(iflag,t13,tmps5)

C-6) t23 element - (w =0 at top)

      iflag=0
      call modecut(iflag,t23,tmps6)


C-Now update elements of SGS tensor array
C-t33,t13 & t23 are zero at top.
      do k=2,nz
        do j=1,ny
          do i=1,nx
            if (k.eq.nz) then

              t11(i,j,k) = tmps1(i,j,k) -
     >                     u(i,j,k)*u(i,j,k)
              t22(i,j,k) = tmps2(i,j,k) -
     >                     v(i,j,k)*v(i,j,k)
              t33(i,j,k) = 0.

              t12(i,j,k) = tmps4(i,j,k) -
     >                     u(i,j,k)*v(i,j,k)
              t13(i,j,k) = 0.
              t23(i,j,k) = 0.

            else

              t11(i,j,k) = tmps1(i,j,k) -
     >                     u(i,j,k)*u(i,j,k)
              t22(i,j,k) = tmps2(i,j,k) -
     >                     v(i,j,k)*v(i,j,k)
              t33(i,j,k) = tmps3(i,j,k) -
     >                     w(i,j,k)*w(i,j,k)
 
              t12(i,j,k) = tmps4(i,j,k) -
     >                     u(i,j,k)*v(i,j,k)
              t13(i,j,k) = tmps5(i,j,k) -
     >                     u(i,j,k)*w(i,j,k)
              t23(i,j,k) = tmps6(i,j,k) -
     >                     v(i,j,k)*w(i,j,k)
            endif
          enddo
        enddo
      enddo

C-Calculate mean temperature. This is a little redundant
C-because calculation will be repeated in postprocessor segment
C-of code.
C-Go to Fourier Space
C-THIS IS DONE BECAUSE WE SEEK THE MEAN VALUE OF THE FILTERED
C-TEMPERATURE FIELD !

      call  horfft (temp,-1,nz)

      do k=1,nz
        tempmean(k)=temp(1,1,k)
      enddo

      call horfft(temp,1,nz)


C-Calculate SGS heat fluxes
C-All SGS fluxes are assumed to obey DIRICHLET BC's at top
C-(Under the assumption that the temperature perturbation is negligible
C-at the top surface)

C-1) th1 element - (theta = 0 at top)

      iflag=0
      call modecut(iflag,th1,tmps1)

C-2) th2 element - (theta  = 0 at top)

      iflag=0
      call modecut(iflag,th2,tmps2)

C-3) th3 element - (theta = 0 at top)

      iflag=0
      call modecut(iflag,th3,tmps3)


C-Now update diagonal elements of SGS heat flux vector
C-Remember: All SGS heatfluxes are zero at top
      do k=2,nz
        do j=1,ny
          do i=1,nx
C-Calculate temperature perturbation of filtered field
            ox(i,j,k) = temp(i,j,k) - tempmean(k)
            if (k.eq.nz) then
              th1(i,j,k) = 0.
              th2(i,j,k) = 0.
              th3(i,j,k) = 0.
            else
              th1(i,j,k) = tmps1(i,j,k) -
     >                     u(i,j,k)*ox(i,j,k)
              th2(i,j,k) = tmps2(i,j,k) -
     >                     v(i,j,k)*ox(i,j,k)
              th3(i,j,k) = tmps3(i,j,k) -
     >                     w(i,j,k)*ox(i,j,k)
            endif
          enddo
        enddo
      enddo
      

      return
      end
