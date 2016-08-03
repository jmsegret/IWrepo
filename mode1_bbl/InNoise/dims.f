      module dims

      IMPLICIT NONE
c-list of variables that indicate indices in each spatial direction
C-6/17/02: There's a change made with regard to nz.
C-nz = total number of points in z-direction (global)
C-1/6/03: Given that the penalty method is discontinuous
C-nz = nzloc*nsubd ...
C-nzloc = number of points in a spectral subdomain.
C-nzm = nzloc-1 (Used in Setup of Legendre interpolant functions etc.)
C-For a continuous method:
C-If we have M subdomains: nz = M*(nzloc-1)+1 .
C-nsubd = number of subdomains in vertical direction = M
c-------------------------------------------------------------------
      Integer , parameter ::  nx = 256, ny=1, nz= 300 
      Integer , parameter ::  nzloc = 25, nsubd = 12 
      Integer , parameter ::  nz2=289, nzaux= 68 

c nsdtype...
c	= 1 --> evenly spaced subdomains
c	= 2 --> geometric distribution
c	= 3 --> read from file
 
      Integer , parameter ::   nsdtype = 1
      Integer , parameter ::   ndims = 2
      Integer , parameter ::   nxpp = nx + 2,  nxh=nx/2, nxhp = nxh+1,
     &              nxf = 3*nx/2 + 1,  nyf = 2*ny,  nyh= ny/2,
     &              nxppy = nxpp*ny,   ntot = nxppy*nz,  
     &              ntot2 = nxppy*nz2, nzm = nzloc - 1
      Integer , parameter ::   nxhpy=nxhp*ny      
c     parameter(nz1=49,nzf=2*nz1)
      Integer , parameter ::   nzfl=96, nzh = nzfl/2 
c-------------------------------------------------------------------

      End module dims 
