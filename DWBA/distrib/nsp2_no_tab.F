      INTEGER FUNCTION NSP2(CH)
C     Find last non-space character in an arbitrary length string.
C
      IMPLICIT NONE
      CHARACTER CH*(*)
      INTEGER I
C
      I = LEN(CH)
913   IF (I.GE.1.AND.CH(I:I).EQ.' ') THEN
          I = I-1
          GO TO 913
      ENDIF
      NSP2 = I
      RETURN
      END
