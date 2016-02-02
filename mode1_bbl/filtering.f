C-This file contains all filtering related subroutines
C************************************************************
C-PHYSICAL SPACE FILTERING ROUTINES (ASSOCIATED WITH
C-ESTIMATION MODEL)
C************************************************************

C ***********************************************************
      Subroutine filtering(isgscalc) 
C-Purpose: To implement 3 way filtering of resolved velocity       
C-field. Retains only larger scales (See Dom, Yee & Luo, TCFD)
C ***********************************************************
C-No change in indexing due to implementation
C-of penalty method required  in this routine.

      include 'dim.h'

      include 'comflow.h'
      include 'commeanv.h'
      include 'comfiltertmp.h'
	include 'comsgs.h'
       

c       open(541,file='mode1')
c       open(542,file='mode2')

       call horfft(u,1,nz)
       call horfft(v,1,nz)
       call horfft(w,1,nz)
       call horfft(temp,1,nz)

c       do i = 1,ntot
c          write(541,*) u(i)
c       end do

C-Before filtering velocity field calculate terms
C-of the type <u_i*u_j> appearing in SGS tensor
        if (isgscalc.eq.1) call sgscalc_pre

C-Iflag =1 -> NEUMANN BC AT TOP
C-Iflag =0 -> DIRICHLET BC AT TOP
       iflag=1
       call modecut(iflag,u,tmps1)
       call modecut(iflag,v,tmps2) 
       iflag=0
       call modecut(iflag,w,tmps3) 
       call modecut(iflag,temp,tmps4) 
    
c      do i = nxppy+1,ntot-nxppy
c        u(i) = tmps1(i)
c        v(i) = tmps2(i)
c        w(i) = tmps3(i)
c        temp(i) = tmps4(i)
c      end do

C-Now update variables with filtered values
C-Special care for k=nz with u,v velocities
       do k=2,nz
c      do k=nzm/8,nz     
         do j=1,ny
            do i=1,nx
               if (k.eq.nz) then
                 u(i,j,k) = tmps1(i,j,k)
                 v(i,j,k) = tmps2(i,j,k)
                 w(i,j,k) = w(i,j,k)
                 temp(i,j,k) = temp(i,j,k)
               else
                 u(i,j,k) = tmps1(i,j,k)
                 v(i,j,k) = tmps2(i,j,k)
                 w(i,j,k) = tmps3(i,j,k)
c                temp(i,j,k) = tmps4(i,j,k)
               endif
            enddo
          enddo
        enddo


C-Now complete SGS calculation
        if (isgscalc.eq.1) call sgscalc_post

c       do i = 1,ntot
c          write(542,*) u(i)
c       end do

       call horfft(u,-1,nz)
       call horfft(v,-1,nz)
       call horfft(w,-1,nz)
       call horfft(temp,-1,nz)

c
	return


      end





C**********************************************************
      subroutine modecut(iflag,u,u0)
C**********************************************************
     
c-Iflag = Flag which indicates whether I consider
c-top boundary and u & v velocities in filterz3pv
c
c This routine obtain large component of the full velocity field by special
c filtering.
c-NOTE: (modification 1/09/03 by PD). In order not to disrupt
c-the original structure of this subroutine

      include 'dim.h'
      include 'comgridptxyz.h'

      dimension u(nxpp,ny,nz),u0(nxpp,ny,nz)
      dimension tmp2(nxpp,ny,nz)
      dimension g3(nxpp,ny,nz),g3s(nxpp,ny,nz),g3c(nxpp,ny,nz) ! ,g3q(ntot)

c filter once
      call filterx3pv(u,tmp2)      !trapezoidal rule
      call filtery3pv(tmp2,g3)
c     call filterz3pv(iflag,tmp2,g3)   ! g3=f(u)
C     call filterz3pv(iflag,u,g3)
      

c filter twice
      call filterx3pv(g3,tmp2)
      call filtery3pv(tmp2,g3s)
c     call filterz3pv(iflag,tmp2,g3s)  ! g3s=f(f(u))
C      call filterz3pv(iflag,g3,g3s)

c filter three times
      call filterx3pv(g3s,tmp2)
      call filtery3pv(tmp2,g3c)
c     call filterz3pv(iflag,tmp2,g3c)  ! g3c=f(f(f(u)))
C     call filterz3pv(iflag,g3s,g3c)

c filter four times
c      call filterx3pv(g3c,g3q)
c      call filtery3pv(g3q,tmp2)
c      call filterz3pv(tmp2,g3q) ! g3q=f(f(f(f(u))))

c obtain large component from full field
      do i=1,nx  
        do j=1,ny
          do k=1,nz
C-Perform operation: u^0=3*f(u)-3*f(f(u))+f(f(f(u)))
            u0(i,j,k)=3.*g3(i,j,k)-3.*g3s(i,j,k)+g3c(i,j,k)  
          enddo
        enddo
      end do

C**PD: Added modification to implement SGS filtering in
C-vertical direction because physical space filter does not work**
      call filter_legensgsz(u0)
      return
      end




C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine filterx3pv(u,uf)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++

c     Filter quantities in x direction using a tophat filter
c     and (3 points) trapezoidal rule for constant grid space

c     Input:  u
c     Output: uf

      include 'dim.h'

      parameter(nxp=nx+1)
      dimension u(nxpp,ny,nz),uf(nxpp,ny,nz),up(0:nxp)

      a=1./4.
      b=2./4.
      c=1./4.

c      call zero(uf,ntot)
      do i=1,nxpp
        do j=1,ny
          do k=1,nz 
            uf(i,j,k) = u(i,j,k)
          enddo
        enddo
      end do

      do 100 ks=1,nsubd

        do kloc=1,nzloc
          k = (ks-1)*nzloc + kloc

          do j=1,ny

            do i=1,nx
              up(i)=u(i,j,k)
            end do
            up(0) =u(nx,j,k)
            up(nxp)=u(1,j,k)
c
            do i=1,nx
              uf(i,j,k)=a*up(i-1)+b*up(i)+c*up(i+1)
            end do
c
          end do
        end do

  100 continue

      return
      end





C +add++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine filtery3pv(u,uf)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Filter quantities in the y direction using a tophat filter
c     and (3 points) trapezoidal rule for constant grid space.

c     Input: u
c     Output: uf

      include 'dim.h'

      parameter(nyp=ny+1)
      dimension u(nxpp,ny,nz),uf(nxpp,ny,nz),up(0:nyp)

      a=1./4.
      b=2./4.
      c=1./4.

c      call zero(uf,ntot)

      do i=1,nxpp
        do j=1,ny
          do k=1,nz
            uf(i,j,k) = u(i,j,k)
          enddo
        enddo
      end do


      do 100 ks=1,nsubd

        do kloc=1,nzloc
          k = (ks-1)*nzloc + kloc

          do i=1,nx
c
            do j=1,ny
              up(j)=u(i,j,k)
            end do
            up(0) =u(i,ny,k)
            up(nyp)=u(i,1,k)
c
            do j=1,ny
              uf(i,j,k)=a*up(j-1)+b*up(j)+c*up(j+1)
            end do
c
          end do
        end do

 100  continue

      return
      end





C****************************************************************
      subroutine filterz3pv(iflag,u,uf)
C****************************************************************

c     Filter data in the z direction using a tophat filter
c     and trapezodal rule.  
c
c     Input: u
c     Output: uf
c
c(**PD**): Iflag is added to deal with free surface at
c          top boundary.
c-Iflag=1, I'm dealing with u,v velocities and consider top surface.
c-Iflag=0,      -""-   with w,T velocities and do not consider top surface.


      include 'dim.h'

      include 'comfilter.h'

      dimension u(nxpp,ny,nz),uf(nxpp,ny,nz)


c      call zero(uf,ntot)
      do i=1,nxpp
        do j=1,ny
          do k=1,nz
            uf(i,j,k) = u(i,j,k)
          enddo
        enddo
      end do

c-kmax=nzm for solid wall at top, k=nz for free surface at top

c     do k=2,nzm

      do 100 ks=1,nsubd

        do 200 kloc=1,nzloc
          k = (ks-1)*nzloc + kloc
          kg = (ks-1)*(nzloc-1) + kloc

          do j=1,ny

            do i=1,nx

C-When working with u or v velocities, assume that velocity above
C-surface is equal to that of surface (Reconsider this when studying
C-subsurface effects !).


C-Filter points in interior of subdomain
            if ((kloc.gt.1).and.(kloc.lt.nzloc)) then

              uf(i,j,k) = atpz(kg)*u(i,j,k-1) + btpz(kg)*u(i,j,k) +
     >                   ctpz(kg)*u(i,j,k+1) 

            else
C-Filter bottom point of subdomain for all subdomains > 1
              if ((ks.gt.1).and.(kloc.eq.1)) then
               
                uf(i,j,k) = atpz(kg)*u(i,j,k-2) + btpz(kg)*u(i,j,k) + 
     >                   ctpz(kg)*u(i,j,k+1)
              endif

C-Filter top point of subdomain for all subdomains < NSUBD
         
              if ((ks.lt.nsubd).and.(kloc.eq.nzloc)) then 
               
                uf(i,j,k) = atpz(kg)*u(i,j,k-1) + btpz(kg)*u(i,j,k) + 
     >                   ctpz(kg)*u(i,j,k+2)
              endif

C-Special treatment for top point at top subdomain:
C-IFLAG=1: If  working with u or v velocities, assume that velocity above
C-surface is equal to that of surface (Reconsider this when studying
C-subsurface effects !).
C-IFLAG=0: Otherwise leave untouched 
C-NOTE: The bottom point (i.e. bottom physical boundary) of the entire
C-domain is not touched.
  
              if ((ks.eq.nsubd).and.(kloc.eq.nzloc))  then

                if (iflag.eq.1) then
                 uf(i,j,k) = 0.25*u(i,j,k-1) + 0.5*u(i,j,k) +
     >                   0.25*u(i,j,k)
                else
                 uf(i,j,k) = u(i,j,k)
                endif

              endif
              
            endif

            enddo               

          enddo

 200    continue
 100  continue

      return
      end


C***************************************************************
C* S P E C T R A L   F I L T E R I N G   S U B R O U T I N E S *
C***************************************************************




C******************************************************************
      subroutine specfilter_uvt(ybar)
C******************************************************************
C
C-Purpose: Called from main to perform Legendre filtering on
C-u,v & T variables
C-PD-6/26/02: THIS IS USED ONLY THE GLOBAL VERTICAL DISCRETIZATION
C-IN THE SUBDOMAIN DECOMPOSITION METHOD, SPECTRAL FILTERING IS PERFORMED
C-IMMEDIATELY UPON COMPUTATION OF THE LEFT & RIGHT HOMOGENEOUS SOLUTION
C-AND PARTICULAR SOLUTION PRIOR TO PATCHING (See Boyd)
C------------------------------------------------------------------
       include 'dim.h'

       include 'paramvarfilterleg.h'
       include 'comflow.h'

       dimension finu(nz),finv(nz),fint(nz),
     > foutu(nz),foutv(nz),foutt(nz)

       call horfft(u,1,nz)
       call horfft(v,1,nz)
       call horfft(temp,1,nz)

       write(*,*) 'FILTERING u,v & T IN LEGENDRE SPACE'
       do i=1,nxpp
         do j=1,ny
c      i=22
c      j=25
C-Set up input function for column
           do k=1,nz
             finu(k) = u(i,j,k)
             finv(k) = v(i,j,k)
             fint(k) = temp(i,j,k)
           enddo
C-Call filtering subroutine
c           call filter_legen(finu,foutu,ybar,iu)
c           call filter_legen(finv,foutv,ybar,iv)
c           call filter_legen(fint,foutt,ybar,it)
C-Now update u,v & T
           do k=1,nz
             u(i,j,k) = foutu(k)
             v(i,j,k) = foutv(k)
             temp(i,j,k) = foutt(k)
           enddo
         enddo
       enddo

       call horfft(u,-1,nz)
       call horfft(v,-1,nz)
       call horfft(temp,-1,nz)

       return
       end



C******************************************************************
      subroutine filter_legen(fin,fout,ybar,ivar,iprt,ks)
C******************************************************************
C------------------------------------------------------------------
C-Purpose to filter out higher Legendre modes and associated numerical
C-noise. (P.D. 4/23/02)
C-Employs basis recombination technique of Boyd and a recommended
C-Legendre-space filter (use exponential for now, gradually transition
C-to Boyd-Vandeven filter).
C-We want to calculate filtered Legendre coefficients in one loop
C-from filtered coefficients of modified base functions.
C------------------------------------------------------------------
C-fin: Input column of data
C-fout: Output column of data (Filtered)
C-ivar: Flag which indicates what variable we are working with
C------------------------------------------------------------------
C-a: Unfiltered Legendre Coefficients
C-af: Filtered    -""-       -""-
C-b: Unfiltered modified base function Coefficients
C-bf: Filtered     -""-   -""- -""-       -""-
C------------------------------------------------------------------
C-NOTE: When treating a function with inhomogeneous BC's one needs
C-To treat the input function with a corresponding extension function
C-(See Boyd's book, pg. 214).
C-
C-PD-6/26/02. SUBDOMAIN SPECTRAL METHOD: All nz-related indices
C-are modified to accommodate calculation in individual subdomains
C-of resolution nzloc.

      logical cond
      logical condsubd
      logical condtemp
      character*14 fileout
      character*7 fmt1

      include 'dim.h'
      include 'paramvarfilterleg.h'
      include 'comlegendre.h'
      include 'comleng.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comspecfilter.h'

      dimension fin(nzloc),fout(nzloc)
      dimension a(0:nzm),b(0:nzm),af(0:nzm),bf(0:nzm)

c     parameter (kc=int(3.*nzm/4.)+1)
c     parameter (kc=int(2.*nzm/3.)+1)
c     parameter (kc=0.50*nzm)
      parameter (pbv = 2.)

      kc = fackc*nzm
      condsubd = ((ks.eq.nsubd).or.(ks.eq.nsubd))
      condtemp = (iprt.eq.1)


C-Use the following when you want to enforce stronger spectral
C-filtering at top & bottom subdomains.
c     if (condsubd) then
c       kc = 0.
c     else
c       kc = 0.5*nzm
c     endif


c     if (iprt.eq.1) then
c       open(49,file='legendre.dat')
c       write(49,*) 'VARIABLES = "K","Z","FIN","FOUT",
c    >  "a_k","a_k^f","b_k","b_k^f","ERROR"'
c     endif

C*******************************************
C-Treat temperature (zeta+ or zeta -) with extension function*
C*******************************************
      if ((ivar.eq.it).or.(ivar.eq.izetap)) then
        do k1=1,nzloc
          zeta = zpts(k1-1)
          fin(k1) = fin(k1) - extend(zeta,it)
        enddo
      endif

      if (ivar.eq.izetam) then
        do k1=1,nzloc          
          zeta = zpts(k1-1)
          fin(k1) = fin(k1) - extend(zeta,izetam)
        enddo
      endif
          
 
C**********************       
C*Begin Filtering Here*
C**********************

C-Calculate Legendre modes (kth mode at point k1)
C-Formulae are in S. 2.3 of Canuto et al. (See also Boyd's
C-book).
      do k=0,nzm

        sum=0.
        do k1=1,nzloc 
          sum = sum + wg(k1-1)*fin(k1)*alpol(k,k1-1)
        enddo

c       write(*,*) alpol(k,nzm)
        if (k.eq.nzm) then
          a(k) = 0.5*nzm*sum
        else
          a(k) = real(k+0.5)*sum
        endif

      enddo

C-Calculate filtered Legendre modes through
C-modes in recombined basis function set and their
C-filtered counterparts.
C-See notes (Wherever in notes you see N, replace it
C-with nzm here).
C

C**********************************************************
C* Filtering for u & v velocities with D-bot & N-top BC's *
C**********************************************************
C-NOTE: Indices of b(k) & bf(k) span 1 to nzm-1 in this case

      cond = ((ivar.eq.iu).or.(ivar.eq.iv))
      if (cond) then

C-Nulify unused modes
      bf(nzm) = 0.
      bf(0) = 0.
C-Highest mode (k=Nzm-1)
      c1 = -float(nzm*nzm)/float((nzm-1)*(nzm-1))

      b(nzm-1) = a(nzm)*c1                       

c     bf(nzm-1) = b(nzm-1)*specfilter_bv(nzm-1,kc,nzm-1,pbv)      
      bf(nzm-1) = b(nzm-1)*specfilter(nzm-1,kc,nzm-1,p,amp)

      af(nzm) = bf(nzm-1)/c1

C-Second highest mode (k=Nzm-2)
      c1 = -float((nzm-1)*(nzm-1))/float((nzm-2)*(nzm-2))
      c2 = float(2*nzm-1)/float(nzm*nzm)

      b(nzm-2) = ( a(nzm-1) - b(nzm-1)*c2)*c1

c     bf(nzm-2) = b(nzm-2)*specfilter_bv(nzm-2,kc,nzm-1,pbv)
      bf(nzm-2) = b(nzm-2)*specfilter(nzm-2,kc,nzm-1,p,amp)


      af(nzm-1) = bf(nzm-2)/c1 + bf(nzm-1)*c2


C-Begin Loop
      do k=nzm-2,2,-1
        c1 = -float(k*k)/float((k-1)*(k-1))
        c2 = float(2*k+1)/float((k+1)*(k+1))
        b(k-1) = ( a(k) - b(k)*c2 - b(k+1))*c1
 
C       bf(k-1) = b(k-1)*specfilter_bv(k-1,kc,nzm-1,pbv)
        bf(k-1) = b(k-1)*specfilter(k-1,kc,nzm-1,p,amp)

        af(k) = bf(k-1)/c1 + bf(k)*c2 + bf(k+1)
  
      enddo

C-(See Boyd), two lowest Legendre modes remain unaffected
      af(1) = a(1)
      af(0) = a(0)

      endif

C**********************************************************
C* Filtering for temperature with D-bot & D-top BC's *
C* (or functions zetap, zetam)
C**********************************************************
C-NOTE: Indices of b(k) & bf(k) span 0 to nzm-2 in this case

      cond = ((ivar.eq.it).or.(ivar.eq.iw).
     >         or.(ivar.eq.izetap).or.(ivar.eq.izetam))
      if (cond)  then            

C-Nulify unused modes
      bf(nzm) = 0.
      bf(nzm-1) = 0.

C-Highest mode (k=Nzm-2)
      b(nzm-2) = a(nzm)

c     bf(nzm-2) = b(nzm-2)*specfilter_bv(nzm-2,kc,nzm-2,pbv)
      bf(nzm-2) = b(nzm-2)*specfilter(nzm-2,kc,nzm-2,p,amp)

      af(nzm) = bf(nzm-2)

c     if (iprt.eq.1) write(*,'(I3,2x,5(e14.4,2x))') 
c    >nzm,specfilter_bv(nzm-2,kc,nzm-2),
c    >a(nzm),af(nzm),b(nzm-2),bf(nzm-2)

C-Second highest mode (k=Nzm-3)
      b(nzm-3) = a(nzm-1)

c     bf(nzm-3) = b(nzm-3)*specfilter_bv(nzm-3,kc,nzm-2,pbv)
      bf(nzm-3) = b(nzm-3)*specfilter(nzm-3,kc,nzm-2,p,amp)


      af(nzm-1) = bf(nzm-3)

c     if (iprt.eq.1) write(*,'(I3,2x,5(e14.4,2x))') 
c    >nzm-1,specfilter_bv(nzm-2,kc,nzm-2),
c    >a(nzm-1),af(nzm-1),b(nzm-3),bf(nzm-3)


C-Begin Loop
      do k=nzm-2,2,-1
        b(k-2) =  a(k) + b(k)

c       bf(k-2) = b(k-2)*specfilter_bv(k-2,kc,nzm-2,pbv)
        bf(k-2) = b(k-2)*specfilter(k-2,kc,nzm-2,p,amp) 

        af(k) = bf(k-2) - bf(k)

c       if (iprt.eq.1) write(*,'(I3,2x,5(e14.4,2x))') 
c    >k,specfilter_bv(k-2,kc,nzm-2),
c    >a(k),af(k),b(k-2),bf(k-2)
 
      enddo

C-(See Boyd), two lowest Legendre modes remain unaffected
      af(1) = a(1)
      af(0) = a(0)

      endif


   
C-Reconstruct function values at each point k1 using kth polynomial
C-Coefficients are filtered Legendre coefficients.
C-If function has been treated with low-order polynomial 
C-function upon input do appropriate output treatment.
      do k1=1,nzloc
        sum = 0.
        do k=0,nzm
          sum = sum+af(k)*alpol(k,k1-1)
        enddo

        fout(k1) = sum

        if ((ivar.eq.it).or.(ivar.eq.izetap)) then 
          zeta = zpts(k1-1)
          fin(k1) = fin(k1)   + extend(zeta,it)
          fout(k1) = fout(k1)  + extend(zeta,it)
        endif

        if (ivar.eq.izetam) then
          zeta = zpts(k1-1)
          fin(k1) = fin(k1)   + extend(zeta,izetam)
          fout(k1) = fout(k1)  + extend(zeta,izetam)
        endif
      enddo

C-Test: Print out modes
c     if (iprt.eq.1) then 
c       write(*,*) 'dumping out to legendre polynomial file'
c       ks = 1
c       do  kloc=1,nzloc
c          k = (ks-1)*(nzloc-1) + kloc
c          zz = z(k)
c          error = abs((fin(kloc)-fout(kloc))/fin(kloc))
c          if (fin(kloc).eq.0) error = 0.
c          write(49,'(1x,I3,2x,9(E12.6,2x))')
c    >     kloc,zz,fin(kloc),fout(kloc),a(kloc-1)**2.,af(kloc-1)**2.,
c    >     b(kloc-1)**2.,bf(kloc-1)**2.,error
c       enddo
c     endif

      close(49)
C-

      return
      end




C******************************************************************
      subroutine filter_legenpen(fin,fout,ybar,iprt,ks)
C******************************************************************
C------------------------------------------------------------------
C-Purpose to filter out higher Legendre modes and associated numerical
C-noise. (P.D. 1/8/03)
C-This is a modified version of the previous Legendre filtering
C-routine adapted to Jan Hesthaven's numerical penalty technique.
C-A simple exponential filter is used (see the review paper
C-by Gottlieb & Hesthaven) and C0-C1 continuity are not enforced
C-as these are not enforced from within the method itself.
C-We want to calculate filtered Legendre coefficients in one loop
C-from filtered coefficients of modified base functions.
C------------------------------------------------------------------
C-fin: Input column of data
C-fout: Output column of data (Filtered)
C------------------------------------------------------------------
C-a: Unfiltered Legendre Coefficients
C-af: Filtered    -""-       -""-
C-b: Unfiltered modified base function Coefficients
C-bf: Filtered     -""-   -""- -""-       -""-
C------------------------------------------------------------------
C-PD-6/26/02. SUBDOMAIN SPECTRAL METHOD: All nz-related indices
C-are modified to accommodate calculation in individual subdomains
C-of resolution nzloc.

      include 'dim.h'
      include 'comlegendre.h'
      include 'comleng.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'        
      include 'comspecfilter.h'

      dimension fin(nzloc),fout(nzloc)
      dimension a(0:nzm),af(0:nzm)
 
      kc = fackc*nzm
C**********************
C*Begin Filtering Here*
C**********************

C-Calculate Legendre modes (kth mode at point k1)
C-Formulae are in S. 2.3 of Canuto et al. (See also Boyd's
C-book).
      do k=0,nzm

        sum=0.
        do k1=1,nzloc
          sum = sum + wg(k1-1)*fin(k1)*alpol(k,k1-1)
        enddo

c       write(*,*) alpol(k,nzm)
        if (k.eq.nzm) then
          a(k) = 0.5*nzm*sum
        else
          a(k) = real(k+0.5)*sum
        endif

      enddo      

C-Calculate filtered Legendre modes
      do k=0,nzm
        af(k) = a(k)*specfilter(k,kc,nzm,p,amp)
c       if (iprt.eq.1) 
c    >  write(*,*) af(k),a(k),specfilter(k,kc,nzm,p,amp)
c       af(k) = a(k)*specfilter_bv(k,kc,nzm,p)
      enddo

C-Reconstruct function values at each point k1 using kth polynomial
C-Coefficients are filtered Legendre coefficients.
      do k1=1,nzloc
        sum = 0.
        do k=0,nzm
          sum = sum+af(k)*alpol(k,k1-1)
        enddo

        fout(k1) = sum
      enddo

      return
      end

C******************************************************************
      subroutine filter_3dlegenfourier(fin3d,ifourier,ilegen,iout)
C******************************************************************
C-(PD: 5/12/03)
C-Routine that filters a 3-D field which is Fourier expanded
C-in the horizontal (x,y).
C-C0-C1 continuity is not required.
      include 'dim.h'
      include 'comlegendre.h'
      include 'comleng.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comwave.h'
      include 'comspecfilter.h'

      parameter (kcf = 0.)
	parameter (tiny=1.e-30)

      dimension fin3d(nxpp,ny,nz)

c     dimension fin(nzloc),fout(nzloc)
      dimension ar(0:nzm),arf(0:nzm)
      dimension ai(0:nzm),aif(0:nzm)
c     dimension out(nxhp,ny),outf(nxhp,ny),filtf(nxhp,ny)

      kc = fackc*nzm
C-Set up file that dumps out data on Fourier transformed
C-field at z=0.5H
c     open (600,file='fourier2dfilt.dat',status='unknown')
c     write(600,*) 'VARIABLES = "k_x","k_y","u","uf","filter"'
c     write(600,*) 'ZONE F=POINT, I=',nxhp,', J=',ny
C**********************
C*Begin Filtering Here*
C**********************
C-PD: THIS OPERATION COULD DEFINITELY BE
C-ACCELERATED BY BEING RECAST INTO SOME MATRIX-MATRIX
C-PRODUCT ! (DO IN FUTURE !)

C-Max. horiz. wave#:
      akmax = sqrt(xsq(nxhp) + ysq(nyh+1))
C-Scan in y-direction
C-Scan in x-direction
      do 10 j=1,ny
        ii = 0
        do 20 i=1,nxpp,2
          ii = ii+1
          ak = sqrt(xsq(ii) + ysq(j))

C-Scan over subdomains
          do 30 ks=1,nsubd
            
C-Focus on ks-subdomain. 
C-If required, first filter 2-D FOURIER modes
            if (ifourier.eq.1) then
c             write(*,*) 'FOURIER'
              do kloc = 1,nzloc
                k = (ks-1)*nzloc+kloc

C-Calculate magnitude and phase of complex 2-D
C-function in Fourier space
                fr = fin3d(i,j,k)
                fi = fin3d(i+1,j,k)
  
c               if ((ii.eq.2).and.(j.eq.2)) write(*,*) fr,fi
                fmag = sqrt(fr*fr+fi*fi)
                cosc = fr/fmag
                sinc = fi/fmag
C-Careful about fourier modes of zero magnitude
	          if (fmag.ne.0.) then
c	          write(*,*) cosc,fr,fi,fmag
                  if (sinc.ge.0.) phase = acos(cosc)
                  if (sinc.lt.0.) phase = -acos(cosc)
C-Now filter according to corresponding wavenumber magnitude

                  cf = specfilter_fourier2d(ak,akmax,kcf,pf,amp)
                  fmagf = fmag * cf
        
C-Reconstruct filtered Fourier function in cartesian form
                  fr = fmagf*cos(phase)
                  fi = fmagf*sin(phase)

			  else
	        
c	            write(*,*) fr,fi,fmag
	            fr = 0.
	            fi = 0.
	 
	          endif
 
                fin3d(i,j,k) = fr
                fin3d(i+1,j,k) = fi
              enddo

            endif

C-Now that Fourier modes have been filtered, proceed with
C-filtering LEGENDRE modes (need to filter both real
C-and imaginary part of Fourier modes)
C-kg represents what k did in previous loop
C-k is a local index
            if (ilegen.eq.1) then
c            write(*,*) 'LEGENDRE'
             do k=0,nzm

               sumr=0.
               sumi=0.
               do k1=1,nzloc
                 kg = (ks-1)*nzloc+k1
                 finr = fin3d(i,j,kg)
                 fini = fin3d(i+1,j,kg)

                 sumr = sumr + wg(k1-1)*finr*alpol(k,k1-1)
                 sumi = sumi + wg(k1-1)*fini*alpol(k,k1-1)
               enddo

c       write(*,*) alpol(k,nzm)
               if (k.eq.nzm) then
                 ar(k) = 0.5*nzm*sumr
                 ai(k) = 0.5*nzm*sumi
               else
                 ar(k) = real(k+0.5)*sumr 
                 ai(k) = real(k+0.5)*sumi
               endif
             enddo

C-Calculate filtered Legendre modes
             do k=0,nzm
               arf(k) = ar(k)*specfilter(k,kc,nzm,p,amp)
               aif(k) = ai(k)*specfilter(k,kc,nzm,p,amp)
             enddo

C-Reconstruct function values at each point k1 using kth polynomial
C-Coefficients are filtered Legendre coefficients.
             do k1=1,nzloc
               kg = (ks-1)*nzloc+k1
               sumr = 0.
               sumi = 0.

               do k=0,nzm
                 sumr = sumr+arf(k)*alpol(k,k1-1)
                 sumi = sumi+aif(k)*alpol(k,k1-1)
               enddo

               fin3d(i,j,kg) = sumr
               fin3d(i+1,j,kg) = sumi
             enddo
           endif
 30       continue           

c       pause
 20     continue
 10   continue

C-Dump out test file
c     do j=nyh+1,ny
c       do i=1,nxhp
c         write(600,'(2x,5(G15.8,1x))')
c    >    xw(i),yw(j),
c    >    out(i,j),outf(i,j),filtf(i,j)
c       enddo
c     enddo

c     do j=1,nyh
c       do i=1,nxhp
c         write(600,'(2x,5(G15.8,1x))')
c    >    xw(i),yw(j),
c    >    out(i,j),outf(i,j),filtf(i,j)
c       enddo
c     enddo

c     close(600)

      return
      end




C******************************************************************
      subroutine filter_legensgsz(uin)
C******************************************************************
C------------------------------------------------------------------
C-Purpose to filter out higher Legendre modes and associated numerical
C-noise. (P.D. 1/8/03)
C-This is a modified version of the previous Legendre filtering
C-routine adapted to Jan Hesthaven's numerical penalty technique.
C-A simple exponential filter is used (see the review paper
C-by Gottlieb & Hesthaven) and C0-C1 continuity are not enforced
C-as these are not enforced from within the method itself.
C-We want to calculate filtered Legendre coefficients in one loop
C-from filtered coefficients of modified base functions.
C-(**PD: 4/1/03**): This specific filter is designed to have a lowerr
C-order and is used as an SGS filter 
C------------------------------------------------------------------
C-u0: 3-D input array
C-fin: Input column of data
C-fout: Output column of data (Filtered)
C------------------------------------------------------------------
C-a: Unfiltered Legendre Coefficients
C-af: Filtered    -""-       -""-
C-b: Unfiltered modified base function Coefficients
C-bf: Filtered     -""-   -""- -""-       -""-
C------------------------------------------------------------------
C-PD-6/26/02. SUBDOMAIN SPECTRAL METHOD: All nz-related indices
C-are modified to accommodate calculation in individual subdomains
C-of resolution nzloc.

      character*14 fileout
      character*7 fmt1

      include 'dim.h'
 
      parameter (psgs=2.)

      include 'comlegendre.h'
      include 'comleng.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comspecfilter.h'
      
      dimension uin(nxpp,ny,nz)
      dimension fin(nzloc) ! ,fout(nzloc)
      dimension a(0:nzm),af(0:nzm)

      kc = fackcbv*nzm


C-Scan over all of physical space and read in data from 3-D array
C-into temporary 1-D array for each subdomain
      do 100 i=1,nx
        do 200 j=1,ny

          do 300 ks=1,nsubd

C-Focus on specific subdomain now and filter within that only

            do kloc=1,nzloc 
              k = (ks-1)*nzloc+kloc
 
              fin(kloc) = uin(i,j,k)
 
            enddo


C**********************
C*Begin Filtering Here*
C**********************
      

C-Calculate Legendre modes (kth mode at point k1)
C-Formulae are in S. 2.3 of Canuto et al. (See also Boyd's
C-book).
            do k=0,nzm

              sum=0.
              do k1=1,nzloc
                sum = sum + wg(k1-1)*fin(k1)*alpol(k,k1-1)
              enddo

c       write(*,*) alpol(k,nzm)
              if (k.eq.nzm) then
                a(k) = 0.5*nzm*sum
              else
                a(k) = real(k+0.5)*sum
              endif


            enddo   

C-Calculate filtered Legendre modes
            do k=0,nzm
c             af(k) = a(k)*specfilter(k,kc,nzm,psgs,amp)
c       if (iprt.eq.1)
c    >  write(*,*) af(k),a(k),specfilter(k,kc,nzm,p,amp)
        af(k) = a(k)*specfilter_bv(k,kc,nzm,psgs)
            enddo

C-Reconstruct function values at each point k1 using kth polynomial
C-Coefficients are filtered Legendre coefficients.
            do k1=1,nzloc
              kgl = (ks-1)*nzloc+k1
              sum = 0.
              do k=0,nzm
                sum = sum+af(k)*alpol(k,k1-1)
              enddo

              uin(i,j,kgl) = sum
c             fout(k1) = sum
            enddo

 300      continue
 200    continue
 100  continue

      return
      end



C******************************************************************
      subroutine speclegencalc(filein,t,fin)
C******************************************************************
C-Diagnostic routine designed to calculate the average Legendre
C-spectra of velocities along energetic vertical centerplane
C-in the center of the flow (j=ny/2).
C-fin: 3-D input function (u,v,w or T)
C-faux: 
      character*1 filein
      character*30 fileout
      character*10 fmt1

      include 'dim.h'
      include 'comlegendre.h'
      include 'comleng.h'
      include 'comspecfilter.h'

      dimension fin(nxpp,ny,nz),faux(nzloc)
      dimension specout(nsubd,nzloc)
      dimension a(0:nzm)


C-Set up output file
      fileout='speclegen'   
      write(unit=fileout(10:10),fmt='(a1)') filein
      fmt1='(f5.3,a4)'
      if (nint(t).ge.10) fmt1 = '(f6.3,a4)'
      if (nint(t).ge.100) fmt1 = '(f7.3,a4)'
      if (nint(t).ge.1000) fmt1 = '(f8.3,a4)'
      write(unit=fileout(11:23),fmt=fmt1) t,'.dat'
      write(*,*) fileout
      open(405,file=fileout)

	


C-Nulify spectral array
      do ks=1,nsubd
        do kloc=1,nzloc
          specout(ks,kloc) = 0.
        enddo
      enddo

C-Convert array back into Physical Space
      call  horfft (fin,1,nz)

	j = nyh

C-Scan over all subdomains
C-(scan layers of subdomains)
      do 100 ks=1,nsubd

        do 200 i=1,nx


          do kloc=1,nzloc
            k = (ks-1)*nzloc + kloc
            faux(kloc) = fin(i,j,k)
c           write(*,*) kloc,faux(kloc)
          enddo  

C-Calculate Legendre modes (kth mode at point k1)
C-Formulae are in S. 2.3 of Canuto et al. (See also Boyd's
C-book).
          do k=0,nzm

            sum=0.
            do k1=1,nzloc
              sum = sum + wg(k1-1)*faux(k1)*alpol(k,k1-1)
            enddo

c       write(*,*) alpol(k,nzm)
            if (k.eq.nzm) then
              a(k) = 0.5*nzm*sum
            else

              a(k) = real(k+0.5)*sum
            endif

          enddo   
c
C-Now calculate Spectra
          do kloc=1,nzloc
            specout(ks,kloc) = specout(ks,kloc) +  
     >                         abs(a(kloc-1))
c           write(*,*) i,ks,kloc,a(kloc-1)
          enddo

 200    continue

C-Average out spectral content
        do kloc=1,nzloc
          specout(ks,kloc) = specout(ks,kloc)/nx
        enddo

 100  continue        

C-Now output to data file 
      do kloc=1,nzloc
        write(405,150) kloc,(specout(ks,kloc),ks=1,nsubd)
      enddo

      close(405)

C-Convert array back into spectral space
      call  horfft (fin,-1,nz)

C-Note: The natuer of this format statement changes with
C-the number of elements in the vertical.
 150  format (1x,I2,2x,10(e12.6,2x))  
 
      return
      end

        
C***********************************************************************
      function specfilter(k,kc,n,p,amp)
C***********************************************************************
C-(PD:4/23/02). Exponential filter in Legendre space.
C-Arguments:
C-a) k = Legendre wavenumber under consideration
C-b) kc = limit wavenumber above which filtering is allowable.
C-c) N = number of points in Legendre expansion

      if (k.ge.kc) then     
       theta = float(k-kc)/float(n-kc)
       specfilter = exp(-amp*theta**(2.*p))
c      specfilter = 1.
       
      else
        specfilter = 1.
      endif

      return
      end



C***********************************************************************
      function specfilter_fourier2d(ak,akmax,kc,p,amp)
C***********************************************************************
C-(PD:4/23/02). Exponential filter in Fourier wave# space.
C-Arguments:
C-a) ak = Fourier wavenumber under consideration
C-b) akmax = Maximum Fourier wavenumber in 2-D Fourier expansion
C-c) kc = fraction of akmax above which filter is applied
C-d) p = filter order
C-e) amp = filter amplitude
      akc = akmax*float(kc)
      
      if (ak.ge.akc) then
       theta = (ak-akc)/(akmax-akc)
c      theta = 0.
       specfilter_fourier2d = exp(-amp*theta**(2.*p))
      else
        specfilter = 1.
      endif

      return
      end


C***********************************************************************
      function specfilter_bv(k,kc,n,p)
C***********************************************************************
C-PD:5/29/02.
C-Implementation of Boyd-Vandeven spectral filter proposed by Boyd (1997)
C-and used by Levin et al. (JCP 1997) in spectral element ocean circulation model.
C-Arguments:
C-a) k = Legendre wavenumber under consideration
C-b) kc = limit wavenumber above which filtering is allowable.
C-This is the lag. Choose kc=2/3*N as indicated in Taylor et al. (JCP 1997)
C-c) N = number of points in Legendre expansion
C-Filter order is set to p=2 to damp out noise at small scales.


      if ((k.gt.kc).and.(k.lt.n)) then

       theta = float(k-kc)/float(n-kc)
c      write(*,*) k,kc,n,theta
C-Proper calculation of chi parameter.
       omega = abs(theta) - 0.5
       if (theta.eq.0.5) then
         chi = 1.
       else
         c1 = 4.*omega**2.
	   c2 = -log(1.-c1)/c1
	   chi = sqrt(c2)
       endif
       arg = 2.*sqrt(p)*chi*omega
       specfilter_bv = 0.5*erfc(arg)

      else
	  if (k.eq.n) then
	    specfilter_bv = 0.
	  else
          specfilter_bv = 1.
	  endif
      endif

      return
      end

C************************************************************************
      function extend(zeta,ief)
C************************************************************************
C-(PD:4/25/02). Extension function that converts function with non-homogeneous
C-BC's to function with homogeneous BC's by simple variable transformation
C-with low-order polynomial.
C-zeta ranges within [-1,1].
C-This function is set up to deal with a pair of D-D BC's.
      include 'paramvarfilterleg.h'
 
      if ((ief.eq.it).or.(ief.eq.izetap)) then
        extend = 0.5 + 0.5*zeta
      endif

      if (ief.eq.izetam)
     > extend = 0.5 - 0.5*zeta
 
      return
      end



C*********************************************************************
      subroutine relaminarize
C*********************************************************************
C-This routine relaminarizes the top and bottom elements to possibly 
C-suppress any numerical
C-instability.
      include 'dim.h' 

      include 'comflow.h'
      include 'commeanv.h' 
      include 'comsubd.h'
      include 'comgridptxyz.h'

      dimension k1(2)
      k1(1) = 1
      k1(2) = nsubd

      do 50 m=1,2
      ks=k1(m)

      write(*,*) 'Relaminarizing element #',m
      do i=1,nx
        do j=1,ny

          do  kloc=1,nzloc
             k = (ks-1)*(nzloc-1) + kloc

C-Set velocity and temperature equal to mean at each point in the
C-top element
             u(i,j,k) = umean(j,k)
             v(i,j,k) = vmean(j,k)
             w(i,j,k) = wmean(j,k)
             temp(i,j,k) = tempmean(k)

c         if ((i.eq.(nxh+1)).and.(j.eq.nyh)) 
c    >    write(*,100) z(k),umean(j,k),vmean(j,k),wmean(j,k),tempmean(k)
c    >    write(*,100) z(k),u(i,j,k),v(i,j,k),w(i,j,k),temp(i,j,k)
         

          enddo  
            
        enddo
      enddo
  50  continue

 100    format (1x,f12.6,2x,4(e12.6,2x))

      return
      end

      


