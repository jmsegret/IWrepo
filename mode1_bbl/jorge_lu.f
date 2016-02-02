!-apv( , ) = is the matrix
!-fvr( ) = the RHS which is the solution when exiting the 2nd routine
!-kpvt( ) = a pivot array. Don't worry about it.
!-nz = is the dimension of the matrix

      dimension apv(nz,nz),kpvt(nz),
     >          fvr(nz)

      call sgetrf(nz,nz,apv,nz,kpvt,info)

      call sgetrs('n',nz,1,apv(1,1),nz,kpvt(1),fvr,nz,info)
