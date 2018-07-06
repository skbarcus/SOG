	SUBROUTINE CHGDE1(A1,A2,A3,A4)
C	This is a routine for the evaluation of the sum of Gaussians (SOG)
C	density defined as in DEPHEL (see for example, PRL 38(77)).
C	The first parameter is a normalization.  Then comes the Gaussian
C	width parameter, then the radii and finally the Gaussian amplitudes
C	For an N-Gaussian SOG, there are therefore 2(N+1) parameters.  The
C	model is taken to be the number of Gaussians.
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION PARAM(48)
	CHARACTER RHONA*(*),RHONAM*26
	COMMON /CHARGE/ NSOG,Z,RP,PARAM
	DATA PI/3.1415926536D0/
	A2=0.D0
	A3=0.D0
	A4=0.D0
C	The following constant is 1/GAMMA**2 in Sick's paper.
	A=3.D0/2.D0/(RP**2)
	B=Z*2.D0*PI/(DSQRT(PI/A)**3)
	DO 205 I=1,NSOG
	    ARG=-PARAM(I)**2*A
	    IF (ARG.GT.-30) THEN
	        B=PARAM(NSOG+I)/(1.D0+2.D0*PARAM(I)**2*A)
	        A2=A2+B*DEXP(ARG)
	        A4=A4+B*PARAM(I)**2*DEXP(ARG)
	    ENDIF
205	CONTINUE
	A4=(8.D0*A**2*A4-2.D0*A*A2)*B/A1**5
	A2=A2*2.D0*B/A1**3
	RETURN
C
	ENTRY CHGDE2(A1,A2,A3)
	R=A2/A1
	A3=0.D0
C	The following constant is 1/GAMMA**2 in Sick's paper.
	A=3.D0/2.D0/(RP**2)
	B=Z*2.D0*PI/(DSQRT(PI/A)**3)
	DO 1 I=1,NSOG
	    C=PARAM(NSOG+I)/(1.D0+2.D0*PARAM(I)**2*A)
	    ARG=-(R-PARAM(I))**2*A
	    IF (ARG.GT.-30) A3=A3+C*DEXP(ARG)
	    ARG=-(R+PARAM(I))**2*A
	    IF (ARG.GT.-30) A3=A3+C*DEXP(ARG)
1	CONTINUE
	A3=A3*B/(A1**3)
	RETURN
C
	ENTRY CHGDE3(RHONA,NPAR)
	NPAR=2*(NSOG+1)
	WRITE (RHONAM,301) NSOG
301	FORMAT ('SOG density (',I2,' Gaussians)')
	RHONA=RHONAM
	RETURN
	END
