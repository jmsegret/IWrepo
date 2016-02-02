C-A common block with parameters relevant to the soliton
C-resuspension problem
C-Essentially **md..(...) is a quantity associated with
C-the wave, not the mean field
       common /soliton/  zeta0,solks,fcshear,alpha0,
     >nmode,isignc
       common /eigen/ eigvf(nz),deigvf(nz),d2eigvf(nz)
       common /eigenuse/fsol(nz),dfsoldz(nz),d2fsoldz2(nz)
       common /sol2dvar/ umwv(nxpp,nz),wmwv(nxpp,nz),
     >                   dumdx(nxpp,nz),dumdz(nxpp,nz),
     >                   dwmdx(nxpp,nz),dwmdz(nxpp,nz),
     >                   rhomwv(nxpp,nz),
     >                   drhomdx(nxpp,nz),drhomdz(nxpp,nz)
