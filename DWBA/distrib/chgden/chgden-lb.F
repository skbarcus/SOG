	SUBROUTINE CHGDE1(A1,A2,A3,A4)
C	This is a charge density routine for a LaGuerre polynomial
C	expansion of the charge distribution.  The first coefficient
C	of the expansion is assumed to be 1 and there is an overall 
C	normalization which can	be calculated by the calling program
C	to force the integral to be Z.  If this feature is not necessary,
C	the normalization may be set to 1.  The parameters are:  1) the
C	normalization, 2) the length parameter, 3) ...  the coefficients.
C	The MODEL number specifies how many coefficients are to be used.
C
C	Converted from the Cyber to Vax 3/6/87
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION COEF(48)
	CHARACTER RHONA*(*),RHONAM*46
	COMMON /CHARGE/ MODEL,RHON,RA,COEF
C
	A2=1.D0
	A3=0.D0
	A4=1.D0
	GL0=1.D0
	NCOEF=MODEL-1
	DO 101 I=1,NCOEF
	    GL0=GL0*(1.D0+.5D0/I)
	    A2=A2+GL0*COEF(I)
	    A4=A4+GL0*COEF(I)*(1.D0+1.33333333333333D0*I)
101	CONTINUE
	A2=A2/A1**3*RHON
	A4=-2.D0*A4/(A1**5*RA**2)*RHON
	RETURN
C
	ENTRY CHGDE2(A1,A2,A3)
	R2=(A2/(A1*RA))**2
	A3=0.D0
	IF (R2.GT.670.D0) RETURN
	GLNM1=1.D0
	GLN=1.5D0-2.D0*R2
	A3=1.D0+COEF(1)*GLN
	IF (NCOEF.GE.2) THEN
	    DO 201 I=2,NCOEF
	        GLNP1=((2*(I-R2)-.5D0)*GLN-(I-.5D0)*GLNM1)/I
	        A3=A3+COEF(I)*GLNP1
	        GLNM1=GLN
	        GLN=GLNP1
201	    CONTINUE
	ENDIF
	A3=A3*EXP(-R2)/A1**3*RHON
	RETURN
C
	ENTRY CHGDE3(RHONA,NPAR)
	NPAR=MODEL+1
	WRITE (RHONAM,301) MODEL
301	FORMAT ('LaGuerre series with overall norm. (',I2,' coef''s)')
	RHONA=RHONAM
	RETURN
	END
