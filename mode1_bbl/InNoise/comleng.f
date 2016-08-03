       module comleng

        

       use dims !, only: nzm
       use comgridptxyz

c      double precision d,d2,d3
       
       IMPLICIT NONE
!      length
!       Real,dimension(*) :: zpts(0:nzm),wg(0:nzm),zpts2(0:64)
!       Real,dimension(*,*) :: d(0:nzm,0:nzm),d2(0:nzm,0:nzm)
!       Real,dimension(*,*) :: d3(0:nzm,0:nzm)
       real :: alp,bet,rv



        


       contains 

!*******************************************************************

       subroutine quad(n,x,w,ndim)
!*******************************************************************
!      include 'dim.h'
!      include 'comlegendre.h'
!      include 'comsubd.h'

!*******************************************************************
! Purpose - to obtain the weighting factors w/wg 

!******************************************************************
! Input
      Integer, INTENT(IN) :: n,ndim
      real, dimension(ndim), INTENT(IN) :: x
   
! output
      real, dimension(ndim), INTENT(OUT) :: w

! parameter
! I am confused as to why we need this 

      Integer, parameter :: nn=1000

! counter 

      INTEGER :: j,k

! Intermediate values 
      real, dimension(nn) :: alp1,al1
!      real, dimension(*) :: xaux(0:nzaux-1),alaux1(0:nn)
      real :: small,xc

!*********************************************************************
c
c  **PD-modify**: Also store values of 0 to ndim Legendre
c  polynomials at all collocation points.
c  determine the Gauss Quadrature weighting factors
!      write(*,*) 'n=', n
c
      small = 1.0e-30
 
      do j = 0 , n
       
          xc = x(j)

!       write(*,*) 'xc = ', xc

       call legen(al1,alp1,n,xc,nn)
c-Calculate weighting factor
       w(j) = 2. /
     +         ( n*(n+1)*al1(n)*al1(n) + small )
c-Store values for j-th Legendre polynomial at point xc=x(k)
c-Store values for derivative of j-th order L. polynomial at
c-collocation point, asll.

!       do j=0,nzm
!         alpol(j,k) = al1(j)
!         dalpol(j,k) = alp1(j)
!       enddo
!        write(*,*) 'w(k) = ', w(j)

      enddo

      end
!************************************************************************


!      write(*,*) 'w=',w 
c     open(210,file='legen_out.dat')
c     do k=1,nz
c       write(210,*) x(k-1),alpol(16,k-1)
c     enddo
c     close(210)

!      goto 100
C-Special addition: Reads in multidomain grid and interpolates
C-values of Legendre polynomials on this grid
C-Read in file with multidomain coordinates
C-USED ONLY IN SINGLE DOMAIN CALCULATIONS
!      open(209,file='z_grid.dat',status='unknown')
!      do k=1,nzaux
!        read(209,*) idummy,xin
!        zaux(k) = xin*.25
!        xaux(k-1) = (2.*xin-zlen)/zlen
c       write(*,*) idummy,xin,xaux(k-1)
!      enddo

C-Now interpolate values of single domain Legendre polynomials
C-PD: After introducing simpler interpolation scheme
C-leave this for future use (4/6/04)
!      do k=0,nzaux-1
c       write(*,*) k,xaux(k)
!        xc = xaux(k)
!        call legenaux(alaux1,n,xc,ndim)

!        do j=0,nzm
!          alpolaux(j,k) = alaux1(j)
!        enddo

!      enddo

c-savva
c     open(210,file='legen_out.dat')
c     do k=1,nzaux
c       write(210,*) xaux(k-1),alpolaux(16,k-1)
c     enddo
c     close(210)

!  100 continue
!      return
      

c**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
c**********************************************************************
C----------------------------------------------------------------------
C-Calculates values of all Legendre polynomials (and immediately highest
C-one in hierarchy) at a given collocation point xc.
C-Same operation for derivatives of Legendre polynomials.

! OUT
!      real, dimension(ndim), INTENT(INOUT) :: al,alp
      Integer, INTENT(IN) :: n,ndim
      real, INTENT(IN) :: xc
      real, dimension(ndim), INTENT(INOUT) :: al,alp

c
! counters 

      integer :: k,km,kp

      al(0) = 1.
      al(1) = xc
      alp(0) = 0.
      alp(1) = 1.
c
      do k=1,n
        kp = k + 1
        km = k - 1
        al(kp) = (2.*k+1.)*xc*al(k)/kp - k*al(km)/kp
      enddo
c
      do k=1,n
        kp = k + 1
        km = k - 1
        alp(kp) = (2.*k+1.)*(al(k)+xc*alp(k))/kp -
     +            k*alp(km)/kp
      enddo
c
      return
      end
!***********************************************************

c**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
c**********************************************************************
c
c  computes the gauss-lobatto collocation points
c
c   N:              degree of approximation (order of polynomials=n+1)
c   ALPHA:          parameter in Jacobi weight
c   BETA:           parameter in Jacobi weight
c   XJAC:           roots from largest to smallest
c
c
c  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
c  for Legendre-Gauss-Lobatto points use           0             0
c
! inputs
      real, INTENT(IN) :: alpha, beta
      INTEGER, INTENT(IN) :: n,ndim
! outputs
      real, dimension(ndim), INTENT(OUT) :: xcol

! intermediate values 
      real, dimension(1000) :: xjac
      INTEGER :: np,i,j,kk,npp,k

      real :: det,rp,rm,a,b,nh,dth,cd,sd,cs,ss,x
      real ::     pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1
      real ::     pnp1m,pdnp1m,pnm,pdnm,pnm1m
      real ::     pnp1, pdnp1, pn, pdn, pnm1

!      real :: rp,rm,a,b,nh,dth,cd,sd,cs,ss,x
      real :: poly,pder,recsum,jm,delx,cssave

!      INTEGER :: i,j
 
! data
      INTEGER :: kstop
      real :: eps 
!      dimension xjac(1000),xcol(0:ndim)
!      common /jacpar/alp,bet,rv
      data kstop/10/
      data eps/1.0e-12/
c

      alp = alpha
      bet = beta
      rv = 1. + alp
      np = n + 1
c

      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1.0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1.0)
      det = pnp*pnm1m - pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m

      a = (rp*pnm1m - rm*pnm1p)/det
      b = (rm*pnp   - rp*pnm)/det
      xjac(1) = 1.0
      nh = (n+1)/2
      dth = 3.14159265/(2*n+1)
      cd = cos(2.*dth)
      sd = sin(2.*dth)
      cs = cos(dth)
      ss = sin(dth)
c
      do 39 j=2,nh
       x = cs
       do 29 k=1,kstop
        call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1 + a*pn + b*pnm1
        pder = pdnp1 + a*pdn + b*pdnm1
        recsum = 0.0
        jm = j - 1
        do 27 i=1,jm
          recsum = recsum + 1.0/(x-xjac(i))
27      continue
28      continue
        delx = -poly/(pder-recsum*poly)
        x = x +delx
        if(abs(delx) .lt. eps) go to 30
29      continue
30      continue
        xjac(j) = x
        cssave = cs*cd - ss*sd
        ss = cs*sd + ss*cd
        cs = cssave
39      continue
        xjac(np) = -1.0
        npp = n + 2
        do 49 i=2,nh
          xjac(npp-i) = -xjac(i)
49      continue
        if(n.ne. 2*(n/2)) then
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         go to 56
        else
        xjac(nh+1) = 0.0
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         return
        endif

c
56      return
        end


c**********************************************************************
        subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
c**********************************************************************

! inputs 
       real, INTENT(IN) :: x        
       INTEGER, INTENT(IN) :: n

! outputs

        real, INTENT(OUT) :: poly,pder,polym1,pderm1,polym2,pderm2

! intermediate values

        real :: apb,polylst,pderlst,a1,a2,a3,a4,b3
        real :: polyn,pdern,psave,pdsave

! counter

        INTEGER :: k

!        write(*,*) 'rv=', rv

!        common /jacpar/ alp,bet,rv
        apb = alp + bet
        poly = 1.0
        pder = 0.0
        if(n.eq.0) return
        polylst = poly
        pderlst = pder
        poly = rv*x
        pder = rv
        if(n.eq.1) return
        do 19 k=2,n
          a1 = 2.*k*(k+apb)*(2.*k+apb-2.)
          a2 = (2.*k+apb-1.)*(alp**2-bet**2)
          b3 = (2.*k+apb-2.)
          a3 = b3*(b3+1.)*(b3+2.)
          a4 = 2.*(K+alp-1.)*(k+bet-1.)*(2.*k+apb)
          polyn = ((a2+a3*x)*poly - a4*polylst) / a1
          pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
          psave = polylst
          pdsave = pderlst
          polylst = poly
          poly = polyn
          pderlst = pder
          pder = pdern
19      continue
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
        return
        end
c














       End module comleng
