	SUBROUTINE CHGDE1(A1,A2,A3,A4)
C	This is a charge density routine for numerically defined
C	densities.  There is only one parameter, the overall
C	normalization. The density is read when the routine is first
C	called using the following format:
C	First card: TITLE -- an 80-character title
C	Then: R0 -- the initial radius in the table of densities
C	      DELTA-R -- the step between entries in the table
C	      R(FINAL) -- the radius of the final entry in the table
C	Then: RHO(R0), RHO(R0+DELTR), RHO(R0+2DELTAR etc. -- The density
C	      table.
C	NOTE: This input is only read when the routine is first called.
C
C	Modified 7/3/81 to run on the NORD-500 (double precision and new
C	calling sequence)
C
C	Modified 8/21/87 for VAX, and to use standard spline interpolation
C	routine.
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	CHARACTER RHONA*(*)
	CHARACTER RHONAM*80
	LOGICAL INIT
	DOUBLE PRECISION R(200),RHO(200),SP(200),B(200),DUMMY(49)
	COMMON /CHARGE/ MODELN,RHON,DUMMY
	DATA INIT/.TRUE./
C
	MTEST=1
	GO TO 1
C
	ENTRY CHGDE2(A1,A2,A3)
	MTEST=2
	GO TO 1
C
	ENTRY CHGDE3(RHONA,NPAR)
	MTEST=3
C
1	IF (INIT) THEN
	    INIT=.FALSE.
	    READ (7,500) RHONAM
500	    FORMAT (A)
	    READ (7,*) R0,DELR,RF
	    NR=(RF-R0)/DELR+1.5D0
	    IF (NR.GT.200) STOP 'Too much input for CHGDEN'
	    READ (7,*) (RHO(I),I=1,NR)
	    DO 900 I=1,NR
	        R(I)=R0+(I-1)*DELR
900	    CONTINUE
	    CALL INSPL2(R,RHO,SP,B,SOM,NR)
	ENDIF
	IF (MTEST-2) 200,300,400
C
200	A2=RHO(1)*RHON/A1**3
	A3=RHON*(RHO(2)-RHO(1))/DELR/A1**4
	A4=(2*RHO(2)-RHO(3)-RHO(1))/DELR**2/A1**5
	RETURN
C
300	RADIUS=A2/A1
	A3=0.D0
	IF (RADIUS.LE.RF) THEN
	    CALL SPLIN2(A2/A1,R,RHO,SP,NR,A3,DRHO,IER)
	    A3=A3*RHON/A1**3
	ENDIF
	RETURN
C
400	RHONA=RHONAM
	NPAR=1
	RETURN
	END
