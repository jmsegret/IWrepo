C-Common block for storage of Legendre polynomials and their derivative
C-functions. First index is polynomial order, second is collocation point
C-index.
C-Also stores points of modified base functions at collocation points.
C-(PD: April 04) Store weights for integration
       common /legendre/alpol(0:nzm,0:nzm),dalpol(0:nzm,0:nzm)
     >,alpol2(0:nzm,0:64),alpolaux(0:nzm,0:nzaux-1)
       common /comlegendre/alpolmod(0:nzm,0:nzm)
 
