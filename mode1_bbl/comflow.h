C-(u,v,w,temp): velocity and temperature fields.
C-(un,vn,wn,tn): contribution to time-stepping scheme
C-at level (n-1).
C-(unm,vnm,wnm,tnm): contribution to time-stepping scheme
C-at level (n-2).
       common  /flow/  u(nxpp,ny,nz),v(nxpp,ny,nz),
     >                 w(nxpp,ny,nz),temp(nxpp,ny,nz),
     >                 un(nxpp,ny,nz),vn(nxpp,ny,nz),
     >                 wn(nxpp,ny,nz),tn(nxpp,ny,nz),
     >                 unm(nxpp,ny,nz),vnm(nxpp,ny,nz),
     >                 wnm(nxpp,ny,nz),tnm(nxpp,ny,nz)
