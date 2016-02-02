C-Contains the two routines needed to perform a rapid
C-LU decomposition of a banded matrix such as that
C-resulting from the penalty solution of the singularly
C-perturbed Helmholtz boundary value problem.


C**********************************************************
C-Created by PD. Stores a band diagonal matrix in compact
C-form and then calls everything necessary to obtain solution.
C-This specific subroutine has been adapted to deal with
C-the band diagonal system encountered in the penalty method.
C-It is set-up for elements with an odd number of points
C-NOTE: This version works with the subroutines from LAPACK
C**********************************************************
      subroutine bandlu_solve(ain,xr,xi,ifuncreal)

      include 'dim.h'
    
      parameter (m=nz,n=nz,kl=nzloc,ku=nzloc)
      parameter (ldab=2*kl+ku+1)
       
      dimension ain(m,n),acomp(ldab,n)
      dimension xr(n),xi(n)
      dimension ipiv(n)
      
C-Store matrix in compact form
C-Follow algorithm recommended in  SGBTRF of LAPACK package
      do j=1,n 
        do i=max(1,j-ku),min(m,j+kl)
      
          acomp(kl+ku+1+i-j,j) = ain(i,j)

        enddo
      enddo
             
C-Perform LU decomposition of compactly stored band matrix
      call sgbtrf(m,n,kl,ku,acomp,ldab,ipiv,info)
      
C-Now solve for real and/or imaginary RHS
C-depending on value of ifuncreal
      if (ifuncreal.eq.1) then
        call sgbtrs( 'n',nz,kl,ku,1,acomp,ldab,ipiv,xr,n,
     >                   info)
      else
        call sgbtrs( 'n',nz,kl,ku,1,acomp,ldab,ipiv,xr,n,
     >                   info)
        call sgbtrs( 'n',nz,kl,ku,1,acomp,ldab,ipiv,xi,n,
     >                   info)
      endif

      return
      end


