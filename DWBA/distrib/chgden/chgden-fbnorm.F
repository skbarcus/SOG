	SUBROUTINE CHGDE1(A1,A2,A3,A4)
C	This is a straight Fourier-Bessel series charge distribution.
C	The COEF's multiply j0 (not like Friar and Negele but like SESFIT)
C	Parameter 1 is the cutoff radius (which should not be fit).
C	Parameters 2, 3, 4,... are coefficients 1, 2, 3, ... The model
C	number should be the number of coefficients to use.  Note that
C	because there is no overall normalization, the coefficients must
C	be chosen to produce a density which integrates to Z.
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION COEF(49)
	CHARACTER RHONA*(*),RHONAM*46
	COMMON /CHARGE/ MODELN,RMAX,COEF
	DATA PI/3.1415926535897D0/
	DATA RMIN/.1D-6/
C
100	A2=1.D0
	A3=0.D0
	A4=-1.D0
	DO 205 I=1,MODELN
	    A2=A2+COEF(I)
	    A4=A4-COEF(I)*DFLOAT(I+1)**2
205	CONTINUE
	A2=A2/A1**3
	RHO0=A2
	A4=A4*(PI/RMAX)**2/A1**5
	RETURN
C
	ENTRY CHGDE2(A1,A2,A3)
	R=A2/A1
	A3=0.D0
	IF (R.GT.RMAX) RETURN
	IF (R.LT.RMIN) GO TO 3
	DO 1 I=1,MODELN
	    ARG=I*PI*R/RMAX
	    A3=A3+COEF(I)*DSIN(ARG)/ARG
1	CONTINUE
	A3=A3/(A1**3)
	RETURN
3	A3=RHO0
	RETURN
C
	ENTRY CHGDE3(RHONA,NPAR)
300	NPAR=MODELN+1
	WRITE (RHONAM,301) MODELN
301	FORMAT ('F.-B. series without overall norm. (',I2,' coef''s)')
	RHONA=RHONAM
	RETURN
	END
