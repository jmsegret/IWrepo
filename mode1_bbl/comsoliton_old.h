C-A common block with parameters relevant to the soliton
C-resuspension problem
C-Essentially **md..(...) is a quantity associated with
C-the wave, not the mean field
       common /soliton/  zeta0,solks,fcshear,nmode,isignc
       common /eigen/ eigvf(nz),deigvf(nz),d2eigvf(nz)
       common /eigenuse/fsol(nz),dfsoldz(nz),d2fsoldz2(nz)
       common /sol2dvar/ umwv(nxpp,nz),wmwv(nxpp,nz),
     >                   dumdx(nxpp,nz),dumdz(nxpp,nz),
     >                   dwmdx(nxpp,nz),dwmdz(nxpp,nz),
     >                   d2umdx2(nxpp,nz),d2umdz2(nxpp,nz),
     >                   d2wmdx2(nxpp,nz),d2wmdz2(nxpp,nz),
     >                   dpmdx(nxpp,nz),dpmdz(nxpp,nz),
     >                   rhomwv(nxpp,nz),
     >                   drhomdx(nxpp,nz),drhomdz(nxpp,nz),
     >                   d2rhomdx2(nxpp,nz),d2rhomdz2(nxpp,nz)
 
