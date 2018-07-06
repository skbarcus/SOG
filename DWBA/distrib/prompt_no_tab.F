      SUBROUTINE PROMPT(IUNIT,STR)
C     Send a character string to the interactive input device leaving
C     the cursor at the end of the line.  This routine nolonger adds
C     extra punctuation at the end of the line (: ? > etc.). That
C     is the responsibility of the calling routine. All trailing
C     spaces are removed and replaced by a single space.
C
      IMPLICIT NONE
      CHARACTER STR*(*)
#     ifdef linux
C         Allocate extra variables when using g77.
          CHARACTER WORK*256
#     endif
      INTEGER IUNIT,NSP2
C
#     ifdef linux
C         g77 does not allow one to concatinate strings with
C         undefined lengths so special work is required.
          WORK = STR(:NSP2(STR))//'"'
          CALL CCOUT(IUNIT,'$',WORK)
#     else
          CALL CCOUT(IUNIT,'$',STR(:NSP2(STR))//'"')
#     endif
      RETURN
      END
