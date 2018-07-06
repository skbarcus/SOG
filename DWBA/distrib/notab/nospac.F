	SUBROUTINE NOSPAC(CH,I1,I2)
C	Find first and last non-space characters in an arbitrary length
C	string.  If string is only spaces, I2 will be 1 less than I1
C
	IMPLICIT NONE
	CHARACTER CH*(*)
	INTEGER I1,NSP1,I2,NSP2
C
	I1 = NSP1(CH)
	I2 = NSP2(CH)
	IF (I2.LT.I1) I2 = I1-1
	RETURN
	END
