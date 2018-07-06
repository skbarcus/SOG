      SUBROUTINE HELLOM(IUNIT,NAME,VER,EDIT,STR)
C     This is a little routine which generates a "hello message"
C     in a standard form which contains the program name (NAME)
C     the version number (VER), the current time and
C     date and the date of last modification (EDIT).  The
C     resulting string is always stored in STR.  If IUNIT
C     is greater than 0, the result is printed on the
C     logical unit IUNIT.
C
      IMPLICIT NONE
      CHARACTER NAME*(*),VER*(*),EDIT*(*),STR*(*)
      CHARACTER TIME*8,DATE*8,OUTPUT*256
      INTEGER IUNIT,NSP1,NSP2
C
      CALL TOD(TIME,DATE)
      WRITE (OUTPUT,10) NAME(NSP1(NAME):NSP2(NAME)),
     + VER(NSP1(VER):NSP2(VER)),TIME,DATE,
     + EDIT(NSP1(EDIT):NSP2(EDIT))
10    FORMAT (A,' version ',A,1X,A,1X,A,' (last edit ',A,')')
      STR = OUTPUT
      IF (IUNIT.NE.0) CALL CCOUT(IUNIT,' ',OUTPUT)
      RETURN
      END
