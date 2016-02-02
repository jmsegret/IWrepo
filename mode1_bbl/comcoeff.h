C-Coefficients to be used in setup of A(n,p)
C-matrix in solution of Helmholtz equation
C-**PD:5/02/02**. Add additional one: diffnbc( , )
C-which deals with Neumann BC's at top surface.
      common /coeff/  diff(nzloc,nzloc), amass(nzloc), 
     >diffnbc(nzloc,nzloc) 

