      Module input


      IMPLICIT None

      Real :: dtin = 0.1, waveperiods = 15.0
      INTEGER :: istart = 0, nstep = 400, tanim = 10000, iform = 0
      Real, parameter :: xlen = 20.0, ylen = 1.0, zlen = 5.0
      Real, parameter :: brunt = 0.2, xnu = 1.e-5, xkappa = 1.e-5
      Real, parameter :: grav = 9.81, expan = 3.38e-3, umax = 0.01
      Real, parameter :: rho0 = 1.e3, xacrit = 2.5, xbcrit = 0.84
      Real, parameter :: zcen =0.13, xacu = 1.7, zsh = 0.2625

      INTEGER :: irelax = 0, ispecfilternl = 1
      INTEGER :: ispecfilterexp = 1, nspecfilterexp = 1
      Real :: amp = 35.0, p = 5.0, pf =10.0, fackc = 0.0
      Real :: fackcbv = 0.0, facdir = 1.e0, facneumn = 1.e0
      Real :: facrobin = 1.e4, tstart = 0.0
      INTEGER :: navgintrfc = 1, ncalclegspec = 10000, ifilter = 0
      INTEGER :: nstepfilter = 5, npostp = 20000, itempreset = 0
      INTEGER :: nrestartdump = 200, n3doutput = 200, i0 = 0
       

      End Module input
