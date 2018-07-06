	INTEGER FUNCTION NSP1(CH)
C	Find first non-space character in an arbitrary length string.
C
	IMPLICIT NONE
	CHARACTER CH*(*)
	INTEGER L,I
C
	L = LEN(CH)
	I = 1
912	IF (I.LT.L.AND.CH(I:I).EQ.' ') THEN
	    I = I+1
	    GO TO 912
	ENDIF
	NSP1 = I
	RETURN
	END
