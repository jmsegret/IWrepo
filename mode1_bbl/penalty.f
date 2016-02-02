C*************************************************************************
C-Contains all subroutines associated with penalty method solution
C-of singularly perturbed 1-D boundary value problem
C-All borrowed from Jan Hesthaven's trilogy in SIAM Journal of Scientific
C-Computing in 1996 and 1997 and his faxed notes.
C*************************************************************************

C************************************************************************
      subroutine penalty_setup(epsilon)
C************************************************************************
C-Calculates whatever penalty method parameters are necessary.
C-These are basically the coefficients associated
C-with the interior subdomains, as the BC's for these subdomains
C-are always Robin (Interior subdomains exist when I have 3 or
C-more subdomains)
      include 'dim.h'

      include 'parambctype.h'
      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comparampen.h'

c     write(*,*) 'Setting up all universal penalty parameters'
C-Set constant omega (All subdomains have the same number of points)
      omega = 2./((nzloc+1)*(nzloc+2))
c     omega = 1./float((nzloc+1)**2)
c     write(*,*) 'omega',omega
C***********************************
C-ALPHA, BETA, GAMMA & DELTA
C-FOR INTERIOR SUBDOMAINS
C***********************************

      do ks=2,nsubd-1
c       write(*,*) 'INTERIOR SUBDOMAINS'
        ybar = zh(ks)
        call param_set(irobin,alpha0,beta0,tau0,epsilon,ybar,ks)
        alpha(ks) = alpha0
        beta(ks) = beta0
        tau1(ks) = tau0

        call param_set(irobin,alpha0,beta0,tau0,epsilon,ybar,ks)
        gamma(ks) = alpha0
        delta(ks) = beta0
        tau2(ks) = tau0

!       write(*,*) 'Subdomain #',ks
!       write(*,*) alpha(ks),beta(ks),gamma(ks),delta(ks),
!    >tau1(ks),tau2(ks)
        
      enddo

      return
      end


C******************************************************************
      subroutine param_set(ifbc,alpha0,beta0,tau0,epsilon,ybar,ks)
C******************************************************************
C-Sets the penalty parameters (alpha,beta,tau) for a given boundary
C-condition (physical boundary or subdomain interface)
      double precision ck,dn,t1,t2,t3
      include 'dim.h'
      include 'parambctype.h'
      include 'comparampen.h'
      include 'compencoeff.h'
C-Penalty parameter coefficients from Jan's papers.
C-DIRICHLET: Min. is 4 to get rid of 1/4 factor in Jan's eqns.
C-According to Jan should be even larger. Must find optimal value
C-through experimentation.
C-NEUMANN: There seems to be no interval of variation for this
C-parameter. Keep fixed at value of 1 
C-recommended (epsilon=1.e-6): facdir =1.e4 and facrobin= 2.5*facdir
C-They should be set to accomodate epsilon. When epsilon is large
C-(Not a singularly perturbed problem, keep them fixed)
C-For small epsilon, one may have to adjust according to observations.
C

C-Dirichlet BC
C-Occurs at top & bottom physical boundaries

      if (ifbc.eq.idirichlet) then

        alpha0 = 1.
        beta0 = 0.

        fac = facdir*(2./ybar)**2.
        tau0 = fac*epsilon/
     >        (alpha0*omega**2.)

      endif

C-Neumann BC
C-Occurs primarily at top physical boundary for existing code
      if (ifbc.eq.ineumann) then

        alpha0 = 0.
        beta0 = 1.

        fac = facneumn*2./ybar
        tau0 = fac/(beta0*omega)

!       write(*,*) ks,tau0
      endif


C-Robin BC
C-Occurs at subdomain interfaces
      if (ifbc.eq.irobin) then

        alpha0 = 1.
        beta0 =  1. ! /epsilon

        ck = alpha0*omega/beta0
        dn = 1./(omega*epsilon*beta0)
        t1 = sqrt(ck**2 + epsilon*ck)
c       write(*,*) facrobin,ck,dn,t1,ybar
c       write(*,*) ks,2.*ck - 2.*t1
        t2 = 2.*ck - 2.*t1

C-Remove this line at some point.
C-Modifies treatment of patching conditions
C-when solvin non-singular problems (i.e. pressure)
        if ((1./epsilon).ge.1.e5) then
          facm = facrobin
        else
          facm = 1. !/sqrt(epsilon)
        endif

        t3 = facm*dn*(epsilon + t2)*
c    >       (2./ybar)**2.
     >         2./ybar

        tau0 = real(t3)
c       tau0 = 70200.
c       bob = (epsilon + 2.*ck - 2.*t1)
c       write(*,*) epsilon,ck,t1,bob
c       write(*,*) tau0

c       ck = 1/(2.*omega)*1./epsilon
c       dn = 1./(omega*ck)
c       t1 = sqrt(1.+ck)
c       tau0 = facrobin*dn*(1+ck-t1)*2./ybar
      
c       write(*,*) omega,epsilon,beta0,dn   
        tau0max = dn*(epsilon + 2.*ck + 2.*t1)*2./ybar
c       tau0max = dn*(1+ck+t1)*2./ybar
        if (tau0.gt.tau0max) then 
          write(*,*) 'Robin BC penalty parameter too high !'
          write(*,*) tau0,tau0max
          pause
        endif

      endif


      return
      end



C*******************************************************************
      subroutine penalty_main(epsilon,ifuncreal,
     >                        bcrbot,bcibot,bcrtop,bcitop,
     >                        ifbcbot,ifbctop,rhsinr,rhsini,
     >                        fr,fi)
C*******************************************************************
    
      logical cbchr,cbchi
C
C-Does some preliminary modifications to set up Helmholtz equation
C-in format indicated by Hesthaven
      include 'dim.h'
      include 'parambctype.h' 
c     parameter (nzmat = nsubd*nzloc)

      include 'comgridptxyz.h'
      include 'comsubd.h'
      include 'comparampen.h'

      dimension fr(nz),fi(nz)
      dimension rhsinr(nz),rhsini(nz)
c     dimension fvr(nzmat),fvi(nzmat)

C******************************
C-MODIFY RHS WHEREVER NECESSARY
C******************************
C-Modify all RHS terms to be multiplied by -epsilon.
      do k=1,nz
        rhsinr(k) = -rhsinr(k)*epsilon
        rhsini(k) = -rhsini(k)*epsilon
      enddo

C**************************
C-ALPHA, BETA, GAMMA, DELTA
C-AT PHYSICAL BOUNDARIES
C**************************
C-NOTE: Everything below assumes homogeneous BC's
C-The set up can handle Dirichlet non-homog. BC's
C-but not Neumann. This must be taken care of in the future !
C-(**PD: 12/18/02**)      

C-Bottom physical boundary:
      ks = 1
      ybar = zh(1)
      call param_set(ifbcbot,alpha0,beta0,tau0,epsilon,ybar,ks)
      alpha(1) = alpha0 
      beta(1) = beta0
      tau1(1) = tau0
C-Weak C1 coupling at top interface of bottom-most element
C-when more than one element used
      if (nsubd.gt.1) then
        call param_set(irobin,alpha0,beta0,tau0,epsilon,ybar,ks) 
        gamma(1) = alpha0
        delta(1) = beta0
        tau2(1) = tau0 
c       write(*,*) alpha0,beta0
      endif

C-Top physical boundary
      ks = nsubd
      ybar = zh(nsubd)
      call param_set(ifbctop,alpha0,beta0,tau0,epsilon,ybar,ks)
      gamma(nsubd) = alpha0
      delta(nsubd) = beta0
      tau2(nsubd) = tau0
C-Weak C1 coupling at bottom interface of topmost element
C-when more than one element used
       if (nsubd.gt.1) then
         call param_set(irobin,alpha0,beta0,tau0,epsilon,ybar,ks) 
         alpha(nsubd) = alpha0
         beta(nsubd) = beta0
         tau1(nsubd) = tau0
c        write(*,*) alpha0,beta0
       endif

!     write(*,*) 'Alpha,Beta,Gamma,Delta',alpha(1),beta(1),
!    >           gamma(nsubd),delta(nsubd)
!     write(*,*) 'Tau1, Tau2',tau1(1),tau2(nsubd)
!     write(*,*) 'Tau1,Tau2 at top/bottom interfaces',
!    >tau2(1),tau1(nsubd)


C************************************
C-NOW CALL HELMHOLTZ EQN. SOLVER
C************************************
       call helmv_penalty(fr,fi,epsilon,rhsinr,rhsini,
     >                    bcrtop,bcitop,bcrbot,bcibot,
     >                    ifbctop,ifbcbot,ifuncreal) 

      return
      end


c********************************************************************
      subroutine helmv_penalty(fvr,fvi,epsilon,rhsinr,rhsini,
     >                         bcrtop,bcitop,bcrbot,bcibot,
     >                         ifbctop,ifbcbot,ifuncreal)
c*********************************************************************
c-Generic Helmholtz equation solver for a vertical column of data.
c-This is a COLLOCATION SOLVER !
c-SPECIALIZED VERSION: Employs Jan Hesthaven's PENALTY method
C-ALSO: Once particular solution is found the procedure is followed
C-to calculate "thickened" BC.
c
c-Peter Diamessis (aka Dr. Dubhead)  Dec. 2002.
c
c---------------------------------------------------------------------
c- Equation solved: f - epsilon * D^2*f = rhsin
c
c
c- epsilon: VERY SMALL COEFFICIENT equal to inverse of
C- k^2 + 2/(nu*dt) or some variant thereof
c-
c- rhsin(z): the initial form of the RHS of the equation at each vertical
c- grid point. Either 0, fixed value
c- or a function which contains contribution from convective terms.
c- NOTE: It has a real and imaginary part. Has also been modified
c- to accomodate Hesthaven's approach.
c-

c-Input: z,,coeffbeta,rhsin,bctop,bcbot,ifbctop,ifbcbot
c
c- bctop = value of boundary condition at top
c- bcbot = value of boundary condition at bottom
c- NOTE: These have a real and imaginary part.
c-
c- ifbcbot, ifbctop = flags which indicate whether the BC at bottom and top,
c- respectively is Dirichlet (Fixed value of f, flag is value 0) or
c- Neumann (Fixed value of dfdz, flag is value 1)
c
c- ifuncreal = Flag that indicates whether I'm working with a real or
c- complex function.
c
c-Output: f
c
c- fv(z) = 1-D vertical array of unknown (solved for) data. Has nz
c- elements. NOTE: It has a real (fvr) and imaginary (fvi) part.
c- Multidomain treatment has motivated transferral to fr & fi


c
c---------------------------------------------------------------------
      include 'dim.h'
      parameter (niter = 500,iprt=0,relerr=10e-4,ido=0)
c     parameter (nzmat = nsubd*nzloc)
      parameter (nmax = 20000)
      parameter (thresh =1.e-14,itol = 4, tol = 1.e-4)
      include 'parambctype.h'

      include 'comcoeff.h'
      include 'comleng.h'
      include 'commatstore.h'
      include 'combcg.h'
      include 'comsubd.h'
      include 'comparampen.h'
      include 'comgridptxyz.h'
      integer ija



      dimension rhsinr(nz),rhsini(nz)

      dimension apv(nz,nz),kpvt(nz),
     >          fvr(nz),fvi(nz)
c     dimension a(nzloc,nzloc),vl(nzloc),vr(nzloc),
c    >          xr(nzloc),xi(nzloc)
c

c
c

c-Preliminary nulifications
       do k5=1,nz
         fvr(k5)=0.
         fvi(k5)=0.
       enddo

C-Nulify linear matrix
C-This shouldn't pose a huge problem
       do k1=1,nz
         do k2=1,nz 
           apv(k1,k2) = 0.
         enddo
       enddo


c      do 100 ks=1,nsubd
c        write(*,*) 'Building matrix for subdomain #',ks
C-Scan all subdomains
C-Set up submatrix and RHS-vector corresponding to specific
C-subdomain. Then plug into main matrix and main RHS vector.
         call matrix_elmsetup(bcrtop,bcrbot,bcitop,bcibot,
     >                           ifbcbot,ifbctop,
     >                           epsilon,rhsinr,rhsini,
     >                           apv,fvr,fvi)

c 100  continue

c      open(900,file='matrix.dat',status='unknown')
c      do i=1,nz
c        write(900,'(15(1x,E11.5))') (apv(i,j),j=1,nz)
c      enddo    
c      close(900)
C-Print patching at bottom row of element
c        do k=1,nz
c          kr2 = nzloc+1 
c          kr3 = 2.*nzloc+1
c          kr4 = 3.*nzloc+1
c          write(*,'(3(2x,E11.5))') apv(kr2,k),apv(kr3,k),apv(kr4,k)
c        enddo
c  Perform LU decomposition of matrix Apv !
c  Use Cray scientific lib.  subroutines
c  First decompose matrix Anp into upper and lower diagonal parts.
c
C-Comment out if using BANDED LU METHOD
c      call sgetrf(nz,nz,apv,nz,kpvt,info)
c      if(info.ne.0) print *,'error in Helmholtz Solver, k= ',info

c

c  Now solve LU decomposed system of equations
c  Use Cray scientific lib.  subroutines
c  Variation if function solved for is complex or real.
c
c
c      if (ifuncreal.eq.1) then
c         call sgetrs('n',nz,1,apv(1,1),nz,kpvt(1),fvr,nz,info)
c       else
c         call sgetrs('n',nz,1,apv(1,1),nz,kpvt(1),fvr,nz,info)
c         call sgetrs('n',nz,1,apv(1,1),nz,kpvt(1),fvi,nz,info)
c      endif

C-Use banded LU method
       call bandlu_solve(apv,fvr,fvi,ifuncreal)

       return
       end




C******************************************************************
      subroutine matrix_elmsetup(bcrtop,bcrbot,bcitop,bcibot,
     >                           ifbcbot,ifbctop,
     >                           epsilon,rhsinr,rhsini,
     >                           a,fvr,fvi)
C******************************************************************
C-Sets up submatrices for linear system of equations in Helmholtz. eqn.
C-Each submatrix corresponds to one element.
C-Variables:
C-a: Submatrix for each element
C-vl: Contribution from element at left (below)
C-vb: Contribution from element at right (above)
C-ks: index of specific element we are working with
C

      character*36 out
      include 'dim.h'

      include 'parambctype.h'

      include 'comcoeff.h'
      include 'comleng.h'
      include 'comsubd.h'
      include 'comparampen.h'

      dimension a(nz,nz),
     >          fvr(nz),fvi(nz)

      dimension rhsinr(nz),rhsini(nz)
c     dimension out(nz)
 
!     open(813,file='takis.dat',status='unknown')
C-Set element lengths
c     write(*,*) 'ybarl - ybarr',ybar,ybarl,ybarr

C-Scan over
      do 100 ks=1,nsubd

c       do k=1,nz
c         bob(k) = 0.
c       enddo

        ybar = zh(ks)
        if (ks.gt.1) ybarl = zh(ks-1)
        if (ks.lt.nsubd) ybarr = zh(ks+1)
C-km is an index which indicates the top-most diagonal
C-element of the specific subdomain submatrix
C-within the greater matrix.
        km = (ks-1)*nzloc
c       write(*,*) 'KM',km
C-Contribution from diffusive operator
        do k1=1,nzloc
          do k2=1,nzloc

               a(km+k1,km+k2) = -epsilon*((2./ybar)**2.)*
     >                           d2(k1-1,k2-1)

          enddo
        enddo

C-Contribution from u-term (appears in diagonal)
      do k1=1,nzloc
        a(km+k1,km+k1) =  a(km+k1,km+k1)  + 1.
      enddo

C-NOTE: This is the most compact approach.
C-IT COULD BE ACCELERATED IN THE FUTURE WITH SPECIALIZED
C-TREATMENT OF THE PHYSICAL BOUNDARIES (i.e. when I have
C-a Neumann condition, all these extra multiplications/additions
C-are not necessary)
C
C-Contribution from boundaries and/or interfaces
C-alpha*u - beta*u,x operator.
C-The same general operation can be applied to
C-all boundaries and interfaces if alpha,beta,gamma,delta & tau1,tau2
C-are properly defined
C-----------------------------
C-Boundary/Interface  Specifics:
C-----------------------------
C-Dirichlet BC: alpha = 1 & beta = 0
C-Neumann BC: alpha = 0 & beta = 1
C-Robin BC: alpha=beta=1 (Only for cases of subdomain interfaces !)
C-BOTTOM (Left) interface/boundary (contribution from alpha & beta)

       a(km+1,km+1) = a(km+1,km+1) + tau1(ks)*alpha(ks) 
    
       cf = -tau1(ks)*beta(ks)*epsilon*2./ybar
       do k=1,nzloc
         a(km+1,km+k) = a(km+1,km+k) + cf*d(0,k-1)
       enddo  

C-TOP (right) interface/boundary (contribution from gamma & delta)

       a(km+nzloc,km+nzloc) = a(km+nzloc,km+nzloc) + 
     >                        tau2(ks)*gamma(ks)

       cf = +tau2(ks)*delta(ks)*epsilon*2./ybar
       do k=1,nzloc
         a(km+nzloc,km+k) = a(km+nzloc,km+k) + 
     >                      cf*d(nzloc-1,k-1)
       enddo

C-Contribution of g(t) & h(t) terms (See Jan Hesthaven's notes)
C-BOTTOM (Left) interface/boundary (contribution from alpha & beta)
C-(See Jan's notes and mine)
C
C-g(t) function
C-If this is the lower physical boundary, g(t) = 0 or some
C-externally assigned value, skip the following part.
       if ((nsubd.gt.1).and.(ks.ne.1)) then
c        write(*,*) 'bottom interface at',ks

C- u(a^i)
         a(km+1,km+1) = a(km+1,km+1) - 0.5*tau1(ks)*alpha(ks)
C- u,x(a^i)
         cf = +0.5*tau1(ks)*beta(ks)*epsilon*2./ybar
         do k=1,nzloc
           a(km+1,km+k) = a(km+1,km+k) + cf*d(0,k-1)
         enddo
C-Nulify matrix elements that originate from two
C-subsequent contributions
         do k=1,nzloc
           a(km+1,km-nzloc+k) = 0.
         enddo
C- u(b^i-1)
         a(km+1,km) = -0.5*tau1(ks)*alpha(ks)
C- u,x(b^i-1)
         cf = +0.5*tau1(ks)*beta(ks)*epsilon*2./ybarl
         do k=1,nzloc
           a(km+1,km-nzloc+k) =   a(km+1,km-nzloc+k) +
     >                            cf*d(nzloc-1,k-1)
c          bob(km-nzloc+k) = bob(km-nzloc+k) + d(nzloc-1,k-1)
         enddo

c        do k=1,nzloc
c          write(*,*) bob(km-nzloc+k)*cf,a(km+1,km-nzloc+k)
c        enddo
c        write(*,*) '------------'

       endif

C-TOP (Right) interface/boundary (contribution from gamma & delta)
C-(See Jan's notes and mine)
C
C-h(t) function
C-If this is the upper physical boundary, h(t) = 0 or some
C-externally assigned value, skip the following part.
       if ((nsubd.gt.1).and.(ks.ne.nsubd)) then
c        write(*,*) 'Top interface at',ks

C- u(b^i)
         a(km+nzloc,km+nzloc) = a(km+nzloc,km+nzloc) - 
     >                    0.5*tau2(ks)*gamma(ks) 
C- u,x(b^i)
         cf = -0.5*tau2(ks)*delta(ks)*epsilon*2./ybar
         do k=1,nzloc
           a(km+nzloc,km+k) = a(km+nzloc,km+k) + cf*d(nzloc-1,k-1)
         enddo
C-Nulify matrix elements that originate from two
C-subsequent contributions
         do k=1,nzloc
           a(km+1,km+nzloc+k) = 0.
         enddo
C- u(a^i+1)
         a(km+nzloc,km+nzloc+1) = a(km+nzloc,km+nzloc+1) -
     >                            0.5*tau2(ks)*gamma(ks)
C- u,x(a^i+1)
         cf = -0.5*tau2(ks)*delta(ks)*epsilon*2./ybarr
c        write(*,*) 'cf',cf
         do k=1,nzloc
c        write(*,*) a(km+nzloc,km+nzloc+k),cf
           a(km+nzloc,km+nzloc+k) = a(km+nzloc,km+nzloc+k)
     >                             + cf*d(0,k-1)
         enddo

       endif

C-Now set up the RHS        
C-Contribution from Non-linear terms etc.
       do kloc=1,nzloc
         k = (ks-1)*nzloc + kloc
         fvr(k) = rhsinr(k)
         fvi(k) = rhsini(k)
       enddo

 100  continue

C-Set up RHS for physical boundaries
C-Bottom
      if (ifbcbot.eq.idirichlet) then
        fvr(1) = rhsinr(1) + tau1(1)*bcrbot
        fvi(1) = rhsini(1) + tau1(1)*bcibot
      else
        fvr(1) = rhsinr(1) - tau1(1)*bcrbot*epsilon
        fvi(1) = rhsini(1) - tau1(1)*bcibot*epsilon
      endif
C-Top
      if (ifbctop.eq.idirichlet) then
        fvr(nz) = rhsinr(nz) + tau2(nsubd)*bcrtop
        fvi(nz) = rhsini(nz) + tau2(nsubd)*bcitop
      else
        fvr(nz) = rhsinr(nz) + tau2(nsubd)*bcrtop*epsilon
        fvi(nz) = rhsini(nz) + tau2(nsubd)*bcitop*epsilon
      endif

!      do i=1,nz
!       do j=1,nz
!         if (a(i,j).ne.0.) write(813,*) i,j,a(j,j)
!       enddo
!     enddo

      
!     do k=1,nz
!       write(*,*) fvr(k),fvi(k)
!     enddo

c     do i=1,nz
c       do j=1,nz
c         if (a(i,j).ne.0) then 
c           out(i)(j:j)='X'
c         else
c           out(i)(j:j)='O'
c         endif
c       enddo
c     enddo

c     do i=1,nz
c     write(*,*) (out(i)(j:j),j=1,nz)
c     enddo
 

      return
      end




C******************************************************************
      subroutine interface_avg(fin)
C******************************************************************
C-When requested averages final function values at subdomain interfaces.
C-This is the 1-D case, where the operation is performed
C-only over a column of data.
      include 'dim.h'
      include 'comgridptxyz.h'
      dimension fin(nz)


C-Scan over all subdomains in vertical
      do 300  ks=2,nsubd
        k = (ks-1)*nzloc + 1 
C-At the subdomain interface points take the average
C-of the the contribution from the two adjacent subdomains.
C-Update first point of interior subdomain
c       write(*,*) z(k),fin(k),fin(k-1)
        fin(k) = 0.5*(fin(k) + fin(k-1))
C-Update last point of subdomain below it
        fin(k-1) = fin(k)
c       write(*,*) z(k),fin(k),fin(k-1)

 300  continue

       return
       end

C******************************************************************
      subroutine interface_avg3dold(fin)
C******************************************************************
C-When requested averages final function values at subdomain interfaces.
C-This is the 3-D case, where the operation is performed
C-over the entire value field of a variable
      include 'dim.h'
      include 'comgridptxyz.h'
      dimension fin(nxpp,ny,nz)


c     open(500,file='testuf.dat',form='unformatted',status='unknown')
C-Scan over all horizontal points
      do 100 i=1,nxpp
        do 200 j=1,ny

C-Scan over all subdomains in vertical
          do 300  ks=2,nsubd
            k = (ks-1)*nzloc + 1
C-At the subdomain interface points take the average
C-of the the contribution from the two adjacent subdomains.
C-Update first point of interior subdomain
c       write(*,*) z(k),fin(k),fin(k-1)
            fin(i,j,k) = 0.5*(fin(i,j,k) + fin(i,j,k-1))
C-Update last point of subdomain below it
            fin(i,j,k-1) = fin(i,j,k)
c       write(*,*) z(k),fin(k),fin(k-1)

 300      continue

c         do ks=1,nsubd
c           do kloc = 1,nzloc
c             k = (ks-1)*nzloc + kloc
c             kg = (ks-1)*(nzloc-1) + kloc
c             write(500) fin(i,j,k)
c           enddo
c         enddo

 200    continue
 100  continue

c     close(500)

       return
       end

C******************************************************************
      subroutine interface_avg3d(fin,ic)
C******************************************************************
C-Perform adaptive averaging at subdomain interfaces that
C-display "jumps"
      include 'dim.h'
      include 'comgridptxyz.h'
      dimension fin(nxpp,ny,nz)
      parameter (ct=0.005)
!AMA 14/October/08 changed ct=0.005 to ct=0.0001 to
! make the solutions smoother


      c1 = 0.5  
      c2 = 0.   
C-Work in Physical Space
      call horfft (fin,1,nz)
C-Scan over all horizontal points
      ic = 0
      do 100 i=1,nxpp
        do 200 j=1,ny

C-Scan over all subdomains in vertical
          do 300  ks=2,nsubd
            k = (ks-1)*nzloc + 1
            km = (ks-2)*nzloc + nzloc
C-At the subdomain interfaces investigate whether the
C-difference of the adjacent subdomain endpoints from
C-the average value is too large (This is essentially
C-equivalent to the criterium of Don et al.)
            fsp = fin(i,j,k)
            fsm = fin(i,j,km)
            favg2 = (fsp+fsm)
            c0 = abs(fsp-fsm)/abs(favg2)

C-If constraint violated, then average !
            if ((c0.gt.ct)) then
C-Follow Fincham's idea: Substitute interfacial points
C-and points above and below with weighted average
            
              ffilter = fin(i,j,km-1)*c1 +
     >                  fin(i,j,km)*c2 +
     >                  fin(i,j,k)*c2+
     >                  fin(i,j,k+1)*c1

              
c             fin(i,j,km-1) = ffilter
              fin(i,j,km) = ffilter
              fin(i,j,k) = ffilter
c             fin(i,j,k+1) = ffilter
              ic = ic + 1
            endif

 300      continue

 200    continue
 100  continue

C-Return to Fourier Space
      call horfft (fin,-1,nz)
c     write(*,*) ic,"points averaged"
       return
       end



C******************************************************************
      subroutine param_output(ifbc,alpha0,beta0,tau0,tau0max,
     >                        epsilon,ybar,ks)
C******************************************************************
C-Sets the penalty parameters (alpha,beta,tau) for a given boundary
C-condition (physical boundary or subdomain interface)
      double precision ck,dn,t1,t2,t3
      include 'dim.h'
      include 'parambctype.h'
      include 'comparampen.h'
      include 'compencoeff.h'
C-Penalty parameter coefficients from Jan's papers.
C-DIRICHLET: Min. is 4 to get rid of 1/4 factor in Jan's eqns.
C-According to Jan should be even larger. Must find optimal value
C-through experimentation.
C-NEUMANN: There seems to be no interval of variation for this
C-parameter. Keep fixed at value of 1
C-recommended (epsilon=1.e-6): facdir =1.e4 and facrobin= 2.5*facdir
C-They should be set to accomodate epsilon. When epsilon is large
C-(Not a singularly perturbed problem, keep them fixed)
C-For small epsilon, one may have to adjust according to observations.
C

C-Dirichlet BC
C-Occurs at top & bottom physical boundaries

      if (ifbc.eq.idirichlet) then

        alpha0 = 1.
        beta0 = 0.

        fac = facdir*(2./ybar)**2.
        tau0 = fac*epsilon/
     >        (alpha0*omega**2.)

      endif

C-Neumann BC
C-Occurs primarily at top physical boundary for existing code
      if (ifbc.eq.ineumann) then

        alpha0 = 0.
        beta0 = 1.

        fac = facneumn*2./ybar
        tau0 = fac/(beta0*omega)

      endif


C-Robin BC
C-Occurs at subdomain interfaces
      if (ifbc.eq.irobin) then

        alpha0 = 1.
        beta0 =  1.

        ck = alpha0*omega/beta0
        dn = 1./(omega*epsilon*beta0)
        t1 = sqrt(ck**2 + epsilon*ck)
c       write(*,*) facrobin,ck,dn,t1,ybar
c       write(*,*) ks,2.*ck - 2.*t1
        t2 = 2.*ck - 2.*t1

C-Remove this line at some point.
C-Modifies treatment of patching conditions
C-when solvin non-singular problems (i.e. pressure)
        if ((1./epsilon).ge.1.e3) then
          facm = facrobin
        else
          facm = 1. !/sqrt(epsilon)
        endif

        t3 = facm*dn*(epsilon + t2)*
     >         2./ybar

        tau0 = real(t3)
c       tau0 = 70200.
c       bob = (epsilon + 2.*ck - 2.*t1)
c       write(*,*) epsilon,ck,t1,bob
c       write(*,*) tau0

c       ck = 1/(2.*omega)*1./epsilon
c       dn = 1./(omega*ck)
c       t1 = sqrt(1.+ck)
c       tau0 = facrobin*dn*(1+ck-t1)*2./ybar

c       write(*,*) omega,epsilon,beta0,dn
        tau0max = dn*(epsilon + 2.*ck + 2.*t1)*2./ybar
c       tau0max = dn*(1+ck+t1)*2./ybar
        if (tau0.gt.tau0max) then
          write(*,*) 'Robin BC penalty parameter too high !'
          write(*,*) tau0,tau0max
          pause
        endif

      endif


      return
      end

C******************************************************************
       subroutine penalty_check(epsilon,ifuncreal,
     >                   bcrbot,bcibot,bcrtop,bcitop,
     >                   ifbcbot,ifbctop,rhsinr,rhsini,
     >                   fvr,fvi)
C******************************************************************
C-Check to see if  linear system is valid, particularly if
C-penalized boundary/patching conditions are satisfied.

      character*36 out
      include 'dim.h'

      include 'parambctype.h'

      include 'comcoeff.h'
      include 'comleng.h'
      include 'comsubd.h'
      include 'comparampen.h'

      dimension a(nz,nz),
     >          fvr(nz),fvi(nz)

      dimension rhsinr(nz),rhsini(nz)
      dimension solv(nz)

C-Set element lengths
c     write(*,*) 'ybarl - ybarr',ybar,ybarl,ybarr

      do k=1,nz
        solv(k) = 0.
      enddo

C-Scan over
      do 100 ks=1,nsubd
        ybar = zh(ks)
        if (ks.gt.1) ybarl = zh(ks-1)
        if (ks.lt.nsubd) ybarr = zh(ks+1)
C-km is an index which indicates the top-most diagonal
C-element of the specific subdomain submatrix
C-within the greater matrix.
        km = (ks-1)*nzloc
c       write(*,*) 'KM',km
C-Contribution from diffusive operator
        do k1=1,nzloc
          do k2=1,nzloc

            solv(km+k1) = solv(km+k1) -
     >      epsilon*((2./ybar)**2.)*d2(k1-1,k2-1)*fvr(km+k2)
          enddo
        enddo

C-Contribution from u-term (appears in diagonal)
      do k1=1,nzloc
        solv(km+k1) = solv(km+k1) + fvr(km+k1) 
      enddo

C-NOTE: This is the most compact approach.
C-IT COULD BE ACCELERATED IN THE FUTURE WITH SPECIALIZED
C-TREATMENT OF THE PHYSICAL BOUNDARIES (i.e. when I have
C-a Neumann condition, all these extra multiplications/additions
C-are not necessary)
C
C-Contribution from boundaries and/or interfaces
C-alpha*u - beta*u,x operator.
C-The same general operation can be applied to
C-all boundaries and interfaces if alpha,beta,gamma,delta & tau1,tau2
C-are properly defined
C-----------------------------
C-Boundary/Interface  Specifics:
C-----------------------------
C-Dirichlet BC: alpha = 1 & beta = 0
C-Neumann BC: alpha = 0 & beta = 1
C-Robin BC: alpha=beta=1 (Only for cases of subdomain interfaces !)
C-BOTTOM (Left) interface/boundary (contribution from alpha & beta)

       solv(km+1) = solv(km+1) + tau1(ks)*alpha(ks)*fvr(km+1)

       cf = -tau1(ks)*beta(ks)*epsilon*2./ybar
       do k=1,nzloc
         solv(km+1) = solv(km+1) + cf*d(0,k-1)*fvr(km+k)
       enddo 

C-TOP (right) interface/boundary (contribution from gamma & delta)

       solv(km+nzloc) = solv(km+nzloc) + 
     >                  tau2(ks)*gamma(ks)*fvr(km+nzloc)

       cf = +tau2(ks)*delta(ks)*epsilon*2./ybar
       do k=1,nzloc
         solv(km+nzloc) = solv(km+nzloc) + cf*d(nzloc-1,k-1)*fvr(km+k)
       enddo

C-Contribution of g(t) & h(t) terms (See Jan Hesthaven's notes)
C-BOTTOM (Left) interface/boundary (contribution from alpha & beta)
C-(See Jan's notes and mine)
C
C-g(t) function
C-If this is the lower physical boundary, g(t) = 0 or some
C-externally assigned value, skip the following part.
       if ((nsubd.gt.1).and.(ks.ne.1)) then
c        write(*,*) 'bottom interface at',ks

C- u(a^i)
         solv(km+1) = solv(km+1) - 0.5*tau1(ks)*alpha(ks)*fvr(km+1)
C- u,x(a^i)
         cf = +0.5*tau1(ks)*beta(ks)*epsilon*2./ybar
         do k=1,nzloc
           solv(km+1) = solv(km+1) +  cf*d(0,k-1)*fvr(km+k)
         enddo

C- u(b^i-1)
         solv(km+1) = solv(km+1) - 0.5*tau1(ks)*alpha(ks)*fvr(km)
C- u,x(b^i-1)
         cf = +0.5*tau1(ks)*beta(ks)*epsilon*2./ybarl
         do k=1,nzloc
           solv(km+1) = solv(km+1) + cf*d(nzloc-1,k-1)*fvr(km-nzloc+k) 
         enddo

       endif

C-TOP (Right) interface/boundary (contribution from gamma & delta)
C-(See Jan's notes and mine)
C
C-h(t) function
C-If this is the upper physical boundary, h(t) = 0 or some
C-externally assigned value, skip the following part.
       if ((nsubd.gt.1).and.(ks.ne.nsubd)) then
c        write(*,*) 'Top interface at',ks

C- u(b^i)
         solv(km+nzloc) = solv(km+nzloc) - 
     >                    0.5*tau2(ks)*gamma(ks)*fvr(km+nzloc)
C- u,x(b^i)
         cf = -0.5*tau2(ks)*delta(ks)*epsilon*2./ybar
         do k=1,nzloc
           solv(km+nzloc) = solv(km+nzloc) + cf*d(nzloc-1,k-1)*fvr(km+k)
         enddo
C- u(a^i+1)
         solv(km+nzloc) = solv(km+nzloc) -
     >                            0.5*tau2(ks)*gamma(ks)*fvr(km+nzloc+1)
C- u,x(a^i+1)
         cf = -0.5*tau2(ks)*delta(ks)*epsilon*2./ybarr
c        write(*,*) 'cf',cf
         do k=1,nzloc
c        write(*,*) a(km+nzloc,km+nzloc+k),cf
           solv(km+nzloc) = solv(km+nzloc)
     >                             + cf*d(0,k-1)*fvr(km+nzloc+k)
         enddo

       endif

 100  continue

      do ks=1,nsubd
        do kloc=1,nzloc
          k =(ks-1)*nzloc+kloc
          write(*,*) solv(k),rhsinr(k)
        enddo
      enddo
          
C-Set up RHS for physical boundaries
C-Bottom
      fvr(1) = rhsinr(1) + tau1(1)*bcrbot
      fvi(1) = rhsini(1) + tau1(1)*bcibot
C-Top
      fvr(nz) = rhsinr(nz) + tau2(nsubd)*bcrtop
      fvi(nz) = rhsini(nz) + tau2(nsubd)*bcitop

c     do i=1,nz
c       do j=1,nz
c         if (a(i,j).ne.0) then
c           out(i)(j:j)='X'
c         else
c           out(i)(j:j)='O'
c         endif
c       enddo
c     enddo

c     do i=1,nz
c     write(*,*) (out(i)(j:j),j=1,nz)
c     enddo


      return
      end


C******************************************************************
      function gmod(zbot,ztop,bcbot,bctop,z)
C******************************************************************
C-Modifies RHS in case of non-homogeneous BC's
    
      zspan = zbot - ztop 
      gmod = (bcbot - bctop)*z/zspan + 
     >       (bctop*zbot - bcbot*ztop)/zspan

      return
      end
