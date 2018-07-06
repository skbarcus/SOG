	SUBROUTINE CHGDE1(A1,A2,A3,A4)
C	Special version of CHGDEN for Larry Cardman and the analysis
C	of the NBS absolute Carbon-12 experiment.
C
C	This version of CHGDEN contains the following models
C
C	     Model Number                  Model
C	          1               parabolic modified Gaussian
C	          2               parabolic modified Fermi
C	          3               harmonic oscillator
C	          4               Stanford C-12 without wiggles
C	          5               Stanford C-12 with wiggles
C	          6               Heisenberg et. al. C-13
C	          7               Stanford O-16
C	          8               Li, Sick, Yeaman SI-28
C	          9               Sinha et. al. (Phys. Rev C 7 (1973)
C	                          1930) CA-40
C	          10              modified harmonic oscillator
C
C	When called at the entry point CHGDE1, there are four arguments.
C	The first argument (A1) is the conversion factor for changing
C	the units of distance used in th.calling program to fermi's.
C	The second argument is filled in with the charge distribution at
C	r=0, the third is the first derivative at r=0 and the fourth is
C	the second derivative of the density at r=0.
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	NMOD is the number of models
	PARAMETER (NMOD=10)
	CHARACTER RHONA*(*),RHONAM(NMOD)*50
	INTEGER NPARS(NMOD)
	DOUBLE PRECISION DUMMY(46)
	COMMON /CHARGE/ MODEL,RHON,C,Z,W,DUMMY
	DATA RHONAM/'PARABOLIC MODIFIED GAUSSIAN',
	1 'PARABOLIC FERMI',
	2 'HARMONIC OSCILLATOR',
	3 'STANFORD CARBON 12 WITHOUT WIGGLES',
	4 'STANFORD CARBON 12 WITH WIGGLES',
	5 'HEISENBERG ET AL CARBON 13',
	6 'STANFORD OXYGEN 16',
	7 'LI, SICK, YEAMAN SILICON 28',
	8 'SINHA ET AL CALCIUM 40',
	9 'MODIFIED HARMONIC OSCILLATOR'/
	DATA NPARS/4,4,4,1,1,1,1,1,1,4/
	IF (MODEL.EQ.1) THEN
C	    Model 1 -- Parabolic modified Gaussian charge distribution
	    AA1=DEXP(-(C/Z)**2)
	    B1=AA1+1.0D0
	    A2=1.0D0/B1
	    A3=0.0D0
	    A4=-2.0D0*AA1/((Z**2)*(B1**2))
	    A4=A4+2.0D0*A2*W/C**2
	ELSEIF (MODEL.EQ.2) THEN
C	    Model 2 -- Parabolic Fermi charge distribution
	    AA1=DEXP(-C/Z)
	    B1=AA1+1.0D0
	    A2=1.0D0/B1
	    A3=-AA1/(Z*(B1**2))
	    A4=2.0D0*W*A2/C**2+A3/Z-2.0D0*A2*A3*AA1/Z
	ELSEIF (MODEL.EQ.3) THEN
C	    Model 3 -- Harmonic well charge distribution
C	    C=K*A
C	    Z=(Z-2)/3,Z between 1 and 8
C	    W=(Z-8)*2/15 for d shell
	    A2=1.0D0
	    A3=0.0D0
	    A4=2.0D0*(W-1.0D0)/C**2
	ELSEIF (MODEL.EQ.4) THEN
C	    Model 4 -- Stanford carbon-12 distribution
C	    No free parameters
	    AA=1.649D0
	    BB=1.705D0
	    AS=0.0019D0
	    BS=4.4D0
	    CS=1.1D0
	    AA2=AA**2
	    B2=BB**2
	    BS2=BS**2
	    CS2=CS**2
	    BC2=BS2/CS2
	    BC2E=DEXP(-BC2)
	    B2D=0.75D0*B2
	    A2=1.0D0+AS*BC2E
	    A3=2.0D0*AS*BS*BC2E/CS2
	    A4=8.0D0/3.0D0/B2-2.0D0/AA2+(4.0D0*AS*BC2/CS2-
	1    2.0D0*AS/CS2)*BC2E
	ELSEIF (MODEL.EQ.5) THEN
C	    Model 5 -- Stanford carbon-12 distribution with wiggles.
C	    No free parameters
	    AA=1.649D0
	    BB=1.705D0
	    AS=0.0019D0
	    BS=4.4D0
	    CS=1.1D0
	    DS=-0.0083D0
	    ES=0.33D0
	    FS=1.275D0
	    GS=2.75D0
	    AA2=AA**2
	    B2=BB**2
	    BS2=BS**2
	    CS2=CS**2
	    ES2=ES**2
	    FS2=FS**2
	    GS2=GS**2
	    BC2=BS2/CS2
	    BC2E=DEXP(-BC2)
	    B2D=0.75D0*B2
	    W0=DS*DEXP(-FS2)*(2.D0+4.D0*ES2*(1.D0-FS2)/GS2)
	    A2=1.0D0+AS*BC2E+W0
	    A3=2.0D0*AS*BS*BC2E/CS2
	    A4=8.0D0/3.0D0/B2-2.0D0/AA2+(4.0D0*AS*BC2/CS2-
	1    2.0D0*AS/CS2)*BC2E+DS*DEXP(-FS2)*(-8.D0*ES2-
	2    2.D0*GS2+8.D0*ES2**2*(4.D0*FS2*(1.D0-FS2)-1.D0)/GS2)
	ELSEIF (MODEL.EQ.6) THEN
C	    Model 6 -- Heisenberg et. al. carbon-13 distribution
	    AA=1.635D0
	    BB=1.594D0
	    AS=0.0012D0
	    BS=4.4D0
	    CS=1.1D0
	    AA2=AA**2
	    B2=BB**2
	    BS2=BS**2
	    CS2=CS**2
	    BC2=BS2/CS2
	    BC2E=DEXP(-BC2)
	    B2D=0.75D0*B2
	    A2=1.0D0+AS*BC2E
	    A3=2.0D0*AS*BS*BC2E/CS2
	    A4=8.0D0/3.0D0/B2-2.0D0/AA2+(4.0D0*AS*BC2/CS2-
	1    2.0D0*AS/CS2)*BC2E
	ELSEIF (MODEL.EQ.7) THEN
C	    Model 7 -- Stanford oxygen-16 charge distribution
C	    No free parameters
	    CC=2.608D0
	    ZZ=.513D0
	    WW=-.051D0
	    GS=2.76D0
	    ES=.35D0
	    DS=-.102D0
	    AA1=DEXP(-CC/ZZ)
	    B1=1.0D0/(AA1+1.0D0)
	    A2=(1.0D0+DS)*B1
	    A3=-AA1*B1**2/ZZ
	    A4=2.0D0*WW*B1/CC**2+A3/ZZ-2.0D0*B1*A3*AA1/ZZ-
	1    (2.0D0*ES**2+GS**2)*DS
	ELSEIF (MODEL.EQ.8) THEN
C	    Model 8 -- Li, Sick, Yeaman silicon-28 distribution
C	    No free parameters
	    CC=3.212D0
	    ZZ=.570D0
	    WW=-.122D0
	    AA1=DEXP(-CC/ZZ)
	    B1=AA1+1.0D0
	    A2=1.0D0/B1
	    A3=-AA1/(ZZ*(B1**2))
	    A4=2.0D0*WW*A2/CC**2+A3/ZZ-2.0D0*A2*A3*AA1/ZZ
	ELSEIF (MODEL.EQ.9) THEN
C	    Model 9 -- Sinha et. al. calcium-40 distribution
C	    No free parameters
	    CC=3.766D0
	    ZZ=.586D0
	    WW=-.161D0
	    Q0=3.14D0
	    P0=.43D0
	    AA=.0814D0
	    AA1=DEXP(-CC/ZZ)
	    B1=1.0D0/(AA1+1.0D0)
	    A2=B1+AA
	    A3=-AA1*B1**2/ZZ
	    A4=2.0D0*WW*B1/CC**2+A3/ZZ-2.0D0*B1*A3*AA1/ZZ-
	1    (P0**2/2.0D0+Q0**2)*AA
	ELSEIF (MODEL.EQ.10) THEN
C	    Model 10 -- modified harmonic oscillator charge distribution.
C	    The third (unitless) parameter should be the nuclear charge.
	    ZP=W
	    IF (W.LE.8.0D0) ZP=(W-2)/3
	    A2=1.0D0
	    A3=0.0D0
	    A4=2.0D0*(ZP/Z**2-1.0D0/C**2)
	ELSE
	    GO TO 9000
	ENDIF
	A2=RHON*A2/A1**3
	A3=RHON*A3/A1**4
	A4=RHON*A4/A1**5
	RETURN
C
	ENTRY CHGDE2(A1,A2,A3)
C	When called at the entry point CHGDE2, the routine returns
C	the value of the charge distribution at r=A2 in the variable
C	argument A3.  The first argument is again the factor which
C	changes distance units to fermies
	R=A2/A1
	A3=0.D0
	IF (MODEL.EQ.1) THEN
C	    Model 1 -- Parabolic modified Gaussian charge distribution
	    ARG=(R**2-C**2)/(Z**2)
	    IF (ARG.LT.50.0D0) THEN
	        A3=1.0D0/(1.0D0+DEXP(ARG))
	        A3=A3*(1.0D0+W*(R**2/C**2))
	    ENDIF
	ELSEIF (MODEL.EQ.2) THEN
C	    Model 2 -- Parabolic Fermi charge distribution
	    ARG=(R-C)/Z
	    IF (ARG.LT.50.0D0) A3=(1.0D0+W*(R/C)**2)/(1.0D0+DEXP(ARG))
	ELSEIF (MODEL.EQ.3) THEN
C	    Model 3 -- Harmonic well charge distribution
C	    C=K*A
C	    Z=(Z-2)/3,Z between 1 and 8
C	    W=(Z-8)*2/15 for d shell
	    ARG=(R/C)**2
	    IF (ARG.LT.87.0D0) A3=(1.0D0+W*ARG)*DEXP(-ARG)
	ELSEIF (MODEL.EQ.4) THEN
C	    Model 4 -- Stanford carbon-12 distribution
C	    No free parameters
	    T2=R**2
	    ARG1=T2/AA2
	    ARG2=T2/B2D
	    ARG3=(R-BS)**2/CS2
	    IF (ARG1.LT.50.0D0) THEN
	        A3=(1.0D0+ARG2)*DEXP(-ARG1)
	        IF (ARG3.LT.50.0D0) A3=A3+AS*DEXP(-ARG3)
	    ENDIF
	ELSEIF (MODEL.EQ.5) THEN
C	    Model 5 -- Stanford carbon-12 distribution with wiggles.
C	    No free parameters
	    R2=R**2
	    A3=(1.D0+R2/B2D)*DEXP(-R2/AA2)+AS*DEXP(-(R-BS)**2/CS2)
	    GT=GS*R
	    IF (GT**2.LE.1D-6) THEN
	        A3=A3+W0
	    ELSE
	        FMET=FS-ES*R
	        FPET=FS+ES*R
	        A3=A3+DS/GT**2*(DEXP(-FMET**2)*(GT*DSIN(GT)-
	1        2.D0*FMET*ES*R*DCOS(GT))+DEXP(-FPET**2)*(GT*DSIN(GT)+
	2        2.D0*FPET*ES*R*DCOS(GT)))
	    ENDIF
	ELSEIF (MODEL.EQ.6) THEN
C	    Model 6 -- Heisenberg et. al. carbon-13 distribution
	    T2=R**2
	    ARG1=T2/AA2
	    ARG2=T2/B2D
	    ARG3=(R-BS)**2/CS2
	    IF (ARG1.LT.50.0D0) THEN
	        A3=(1.0D0+ARG2)*DEXP(-ARG1)
	        IF (ARG3.LT.50.0D0) A3=A3+AS*DEXP(-ARG3)
	    ENDIF
	ELSEIF (MODEL.EQ.7) THEN
C	    Model 7 -- Stanford oxygen-16 charge distribution
C	    No free parameters
	    ARG=(R-CC)/ZZ
	    IF (ARG.LT.50.0D0) THEN
	        A3=(1.0D0+WW*(R/CC)**2)/(1.0D0+DEXP(ARG))
	        IF (GS*R.LT.1D-6) THEN
	            A3=A3+DS*B1
	        ELSE
	            A3=A3+B1*DS*DSIN(GS*R)*DEXP(-(ES*R)**2)/(GS*R)
	        ENDIF
	    ENDIF
	ELSEIF (MODEL.EQ.8) THEN
C	    Model 8 -- Li, Sick, Yeaman silicon-28 distribution
C	    No free parameters
	    ARG=(R-CC)/ZZ
	    IF (ARG.LT.50.0D0) A3=(1.0D0+WW*(R/CC)**2)/(1.0D0+DEXP(ARG))
	ELSEIF (MODEL.EQ.9) THEN
C	    Model 9 -- Sinha et. al. calcium-40 distribution
C	    No free parameters
	    ARG=(R-CC)/ZZ
	    IF (ARG.LT.50.0D0) THEN
	        A3=(1.0D0+WW*(R/CC)**2)/(1.0D0+DEXP(ARG))
	        IF (Q0*R.LE.1D-6) THEN
	            A3=A3+AA
	        ELSE
	            A3=A3+AA*DSIN(Q0*R)*DEXP(-(P0*R)**2/4.0D0)/(Q0*R)
	        ENDIF
	    ENDIF
	ELSEIF (MODEL.EQ.10) THEN
C	    Model 10 -- modified harmonic oscillator charge distribution.
C	    The third (unitless) parameter should be the nuclear charge.
	    ARG=(R/C)**2
	    IF (ARG.LT.87.0D0) A3=(1.0D0+ZP*(R/Z)**2)*DEXP(-ARG)
	ELSE
	    GO TO 9000
	ENDIF
	A3=A3*RHON/A1**3
	RETURN
C
	ENTRY CHGDE3(RHONA,NPAR)
C	When called at the entry CHGDE3, a name for the density is
C	returned in the character variable RHONA and the number of
C	parameters used by the density are returned in NPAR.
	IF (MODEL.LT.1.OR.MODEL.GT.NMOD) GO TO 9000
	RHONA=RHONAM(MODEL)
	NPAR=NPARS(MODEL)
	RETURN
9000	STOP 'Bad Model Number -- ABORT!'
	END
