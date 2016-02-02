C-Contains all subroutines needed to fit a Gaussian profile to
C-the mean y/z velocity profiles at the associated centerplanes.


c -----------------------------------------------------------------------
c
c	LSN
c
c	Test subroutine LSFIT, which calls the infamous CURFIT
c
c ------------------------------------------------------------------------
c
	subroutine curvefit(ampl0,shift0,width0,y,x,nin,
     >                      ampl,shift,width)
	IMPLICIT REAL*8(A-H,O-Z)
        parameter (n=1000)
        dimension x(n),y(n),a(10),yfit(n),sigmaa(10),
     >            sigmay(n),deltaa(10)
c	DIMENSION X(25),Y(25),A(10),YFIT(25),SIGMAA(10),SIGMAY(25),DELTAA(10)


C-READ IN DATA FROM FILE (In actual code this should be just an array
C-from common block)
c	data (x(i),i=1,10)/1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
c	data (y(i),i=1,10)/1.,1.,3.,5.,8.,6.,5.,4.,2.,1./

	npts=n
	nterms=3

c	Set up initial guesses

	a(1)=ampl0      !amplitude
	a(2)=shift0     !shift
	a(3)=width0      !width
	a(4)=0.		!dc offset
	a(5)=0.		!linear
	a(6)=0.		!quadratic

        call lsfit(x,y,a,NPTS,NTERMS,yfit,sigmaa,sigmay,deltaa)

c       do i=1,nin
c         write(*,'(1x,3(E12.6,2x))') x(i),y(i),yfit(i)
c       enddo

c       write(*,*)  ampl0,shift0,width0
c       write(*,*)  (a(i),i=1,3)

c       open(500,file='out.dat',status='unknown')
c       write(500,*) 'VARIABLES = "x","y","yfit"'
c       do i=1,nin
c         write(500,'(1x,3(E12.6,2x))') x(i),y(i),yfit(i)
c       enddo
         
        
        ampl = a(1)
        shift = a(2)
        width = a(3)

	return         
	end

C
C	F U N C T N   I
c	Example Gaussian function mimics PV-WAVE function.  Works well
c	when NPTS is large, and .not.well otherwise.
C

c	FUNCTION FUNCTN(X,I,A)

c	IMPLICIT REAL*8(A-H,O-Z)
c	DIMENSION X(25),A(10)

c	z=(x(i)-A(2))/A(3)
c	y=A(1)*exp(-z**2/2.) + A(4) + A(5)*x(i) + A(6)*x(i)**2

c	FUNCTN=y

c	RETURN
c	END

C
C	F U N C T N   I I
c	Simple Gaussian function, whose parameters are:
c	a1 = amplitude
c	a2 = shift (location in x)
c	a3 = sigma (half width)
C

	FUNCTION FUNCTN(X,I,A)

	IMPLICIT REAL*8(A-H,O-Z)
        parameter (n=1000)
	DIMENSION X(N),A(10)

	z=(x(i)-A(2))/A(3)
	y=A(1)*exp(-z**2/2.)

	FUNCTN=y

	RETURN
	END


************************************************************************


C	         See Bevington 
C	         (Least Square Estimation of Nonlinear Parameter)
C	         LSNFIT performs least squares fit on arbitrary functions.
C	         Program does Taylor expansion to linearize the function
C	         and iterates to minimize CHI**2. The following subroutine
C	         must be supplied.

C
C	FUNCTION FUNCTN(X,I,A)
C	DIMENSION X(1,1), A(1)
C	...
C	...
C	RETURN
C	END

C
C	I is the subscript
C	X(1,1) is the independent variable
C	A is the array of parameters
C

C	---------------------------------------------------------------
C	I N P U T  D A T A 
C
C	NPTS, NTERMS, MODE, MAXTRY, ALIM, FLAMDA, X, Y, SIGMAY,
C	A, DELTAA		
C	---------------------------------------------------------------
C
C	NPTS    number of data points to be fitted in this data set
C	NTERMS  number of parameters to be adjusted  ( < NPTS )
C	MODE    1 -- 1/sigma**2 fit---NOTTT !!-now is Y(i) fit...8/3/92
C		0 -- equal weight fit
C	  	-1 -- 1/Y(i) fit
C	MAXTRY  maximum number of iterations to minimize chisqr
C	        (set to 20 for reasonable functions)
C	ALIM    fractional change between 2 successive chisqr before
C		iteration procedure stops.
C		(set to 0.005 for reasonable functions. 0.01 gives 
C		slightly less accuracy.)
C	FLAMDA  parameter that controls iteration procedure
C		(set to 0.01)
C	X(I)  array of independent variables
C	Y(I)    array of dependent variables	
C	SIGMAY(I)  array of one sigma error on Y(I)
C	A(I)    initial guesses for the array of parameters
C		(initial guesses must be supplied)
C	DELTAA(I)  array of differences for which derivatives of
C		   functions with repect to A(I) are calculated
C
C		   dF   F(A(I)+DELTAA(I)) - F(A(I) - DELTAA(I))

C		   -- = ---------------------------------------

C		   dA                 2*DELTAA(I)
C

C	fits up to 25 data points, 3 independent variables and 10 parameters

	Subroutine lsfit(x,y,a,NPTS,NTERMS,yfit,sigmaa,sigmay,deltaa)
	IMPLICIT REAL*8(A-H,O-Z)
        parameter (nmax=1000)
	DIMENSION X(nmax),Y(nmax),A(10),YFIT(nmax),SIGMAA(10),SIGMAY(nmax),
     >  DELTAA(10)
	DIMENSION B(10)
c
c
 	IN = 0
	MODE=1
	MAXTRY = 20
	ALIM = 0.005
	FLAMD = 0.01

10	IN = IN + 1

	do 52 i=1,nterms
52	DELTAA(i)=.0001

	NFREE = NPTS - NTERM

	DO 40 I=1,NPTS
	YFIT(I) = FUNCTN(X,I,A)
40	continue

	CHISQR = FCHISQ(Y,SIGMAY,NPTS,NFREE,MODE,YFIT)
	CHISQO = 0.
	NOTRYS = 0
	FLAMDA = FLAMD

	DO 91 I=1,NTERMS
91	B(I) = A(I)

90	NOTRYS = NOTRYS + 1

	 IF (FLAMDA.LT.1.E-5) FLAMDA=1.E-5


	CALL CURFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,DELTAA,SIGMAA,
     1  FLAMDA,YFIT,CHISQR)

	IF (CHISQR.LT.0.) GO TO 10
	IF (NOTRYS.GT.MAXTRY) GO TO 100

c	IF(CHISQR .EQ. 0.) GO TO 100
C	Use epsilon test
C	DO 92 I=1,NTERMS
C	IF (ABS(B(I)-A(I))/(TAU+ABS(A(I))).GT.EPS) GO TO 93
C92	CONTINUE
C	GO TO 100

	IF ((ABS(CHISQO-CHISQR)/CHISQR).LT.ALIM) GO TO 100
93	CHISQO = CHISQR

	DO 94 I=1,NTERMS
94	B(I) = A(I)

	GO TO 90

100	continue

200	return

	END

C

C	C U R F I T  (Bevington, Ch.11, pp.232-242)

C

	SUBROUTINE CURFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,DELTAA,
     1  SIGMAA,FLAMDA,YFIT,CHISQR)

        parameter (nmax=1000)
	IMPLICIT REAL*8(A-H,O-Z)
C	COMMON/FIT/NTRYS
C	DOUBLE PRECISION ARRAY
	DIMENSION X(nmax),Y(nmax),SIGMAY(nmax),A(10),DELTAA(10),
     1  YFIT(nmax),SIGMAA(10)
	DIMENSION WEIGHT(nmax),ALPHA(10,10),BETA(10),DERIV(10),
     1  ARRAY(10,10),B(10)

11	NFREE = NPTS - NTERMS
	IF (NFREE) 13,13,20
13	CHISQR = 0
	GO TO 110

20	DO 30 I=1,NPTS
21		IF (MODE) 22,27,29
22		IF (Y(I)) 25,27,23
23		WEIGHT(I) = 1./Y(I)
		GO TO 30
25		WEIGHT(I) = 1./(-Y(I))
		GO TO 30
27		WEIGHT(I) = 1.
		GO TO 30
29		if (y(i)) 24,26,28
28		weight(i)=y(i)
		go to 30
24		weight(i)=-y(i)
		go to 30
26		weight(i)=.1
		go to 30
30	CONTINUE
31	DO 34 J=1,NTERMS
		BETA(J) = 0
		DO 34 K=1,J
34			ALPHA(J,K) = 0
41	DO 50 I=1,NPTS
		CALL FDERIV(X,I,A,DELTAA,NTERMS,DERIV)
		DO 46 J=1,NTERMS
		BETA(J) = BETA(J) + WEIGHT(I)*(Y(I)-FUNCTN(X,I,A))*DERIV(J)
		DO 46 K=1,J
46		ALPHA(J,K) = ALPHA(J,K) + WEIGHT(I)*DERIV(J)*DERIV(K)
50	CONTINUE
51	DO 53 J=1,NTERMS
		DO 53 K=1,J
53		ALPHA(K,J) = ALPHA(J,K)
61	DO 62 I=1,NPTS
62		YFIT(I) = FUNCTN(X,I,A)
63	CHISQ1 = FCHISQ(Y,SIGMAY,NPTS,NFREE,MODE,YFIT)
	NTRYS = 0
70	NTRYS = NTRYS + 1
71	DO 74 J=1,NTERMS
		DO 73 K=1,NTERMS
73		ARRAY(J,K) = ALPHA(J,K)/SQRT(ALPHA(J,J)*ALPHA(K,K))
74	ARRAY(J,J) = 1. + FLAMDA
80	CALL MATINV(ARRAY,NTERMS,DET)
81	DO 84 J=1,NTERMS
		B(J) = A(J)
		DO 84 K=1,NTERMS
84		B(J) = B(J) + BETA(K)*ARRAY(J,K)/SQRT(ALPHA(J,J)*ALPHA(K,K))
91	DO 92 I=1,NPTS
92		YFIT(I) = FUNCTN(X,I,B)
93	CHISQR = FCHISQ(Y,SIGMAY,NPTS,NFREE,MODE,YFIT)
	IF (CHISQ1-CHISQR) 95,101,101
95	FLAMDA = 10.*FLAMDA
	IF (NTRYS.LT.10) GO TO 70
101	DO 103 J=1,NTERMS
		A(J) = B(J)
103		SIGMAA(J) = SQRT(ARRAY(J,J)/ALPHA(J,J))
	FLAMDA = FLAMDA/10.
110	RETURN

	END

C
C	M A T I N V
C

	SUBROUTINE MATINV(ARRAY,NORDER,DET)
	IMPLICIT REAL*8(A-H,O-Z)
C	DOUBLE PRECISION ARRAY,AMAX,SAVE
	DIMENSION ARRAY(10,10),IK(10),JK(10)
10	DET = 1.
11	DO 100 K=1,NORDER
	AMAX = 0.
21	DO 30 I=K,NORDER
	DO 30 J=K,NORDER
23	IF (ABS(AMAX)-ABS(ARRAY(I,J))) 24,24,30
24	AMAX = ARRAY(I,J)
	IK(K) = I
	JK(K) = J
30	CONTINUE
31	IF (AMAX) 41,32,41
32	DET = 0.
	GO TO 140
41	I = IK(K)
	IF (I-K) 21,51,43
43	DO 50 J=1,NORDER
	SAVE = ARRAY(K,J)
	ARRAY(K,J) = ARRAY(I,J)
50	ARRAY(I,J) = -SAVE
51	J = JK(K)
	IF (J-K) 21,61,53
53	DO 60 I=1,NORDER
	SAVE = ARRAY(I,K)
	ARRAY(I,K) = ARRAY(I,J)
60	ARRAY(I,J) = -SAVE
61	DO 70 I=1,NORDER
	IF (I-K) 63,70,63
63	ARRAY(I,K) = -ARRAY(I,K)/AMAX
70	CONTINUE
71	DO 80 I=1,NORDER
	DO 80 J=1,NORDER
	IF (I-K) 74,80,74
74	IF (J-K) 75,80,75
75	ARRAY(I,J) = ARRAY(I,J) + ARRAY(I,K)*ARRAY(K,J)
80	CONTINUE
81	DO 90 J=1,NORDER
	IF (J-K) 83,90,83
83	ARRAY(K,J) = ARRAY(K,J)/AMAX
90	CONTINUE
	ARRAY(K,K) = 1./AMAX
100	DET = DET*AMAX
101	DO 130 L=1,NORDER
	K = NORDER - L + 1
	J = IK(K)
	IF (J-K) 111,111,105
105	DO 110 I=1,NORDER
	SAVE = ARRAY(I,K)
	ARRAY(I,K) = -ARRAY(I,J)
110	ARRAY(I,J) = SAVE
111	I = JK(K)
	IF (I-K) 130,130,113
113	DO 120 J=1,NORDER
	SAVE = ARRAY(K,J)
	ARRAY(K,J) = -ARRAY(I,J)
120	ARRAY(I,J) = SAVE
130	CONTINUE
140 	RETURN

	END

C
C	F D E R I V
C
	SUBROUTINE FDERIV(X,I,A,DELTAA,NTERMS,DERIV)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION X(10,10),A(10),DELTAA(10),DERIV(10)
11	DO 18 J=1,NTERMS
	AJ = A(J)
	DELTA = DELTAA(J)
	A(J) = AJ + DELTA
	YFIT = FUNCTN(X,I,A)
	A(J) = AJ - DELTA
	DERIV(J) = (YFIT - FUNCTN(X,I,A))/(2.*DELTA)
18	A(J) = AJ
	RETURN

	END

C

C	F C H I S Q

C

	FUNCTION FCHISQ(Y,SIGMAY,NPTS,NFREE,MODE,YFIT)

	IMPLICIT REAL*8(A-H,O-Z)
C	DOUBLE PRECISION CHISQ, WEIGHT
        parameter (nmax=1000)
	DIMENSION Y(nmax),SIGMAY(nmax),YFIT(nmax)
11	CHISQ = 0
12	IF (NFREE) 13,13,20
13	FCHISQ = 0.
	GO TO 40
20	DO 30 I=1,NPTS
21	IF (MODE) 22,27,29
22	IF (Y(I)) 25,27,23
23	WEIGHT = 1./Y(I)
	GO TO 30
25	WEIGHT = 1./(-Y(I))
	GO TO 30
27	WEIGHT = 1.
	GO TO 30
29	if (y(i)) 24,26,28
28	weight=y(i)
	go to 30
24	weight=-y(i)
	go to 30
26	weight=.1
	go to 30

30	CHISQ = CHISQ + WEIGHT*(Y(I)-YFIT(I))**2
31	FREE = NFREE
32	FCHISQ = CHISQ/FREE
40	RETURN

	END


