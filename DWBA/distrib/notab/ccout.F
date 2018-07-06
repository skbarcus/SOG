	SUBROUTINE CCOUT(IUNIT,CC,STR)
C	This routine prints a character string (STR) with carriage
C	control determined by the character CC on logical unit IUNIT.
C	Trailing blanks are not printed. However, double quote characters
C	are replaced by "non-removable" spaces.
C
C	This version has been modified to allow for the two different
C	calling sequences of FPUTC for Unix.
C
C	Modified to produce correct output if the Absoft Fortran-77
C	compiler is used on a Unix system.
C
	IMPLICIT NONE
C
#	ifdef unix
	    INTEGER IERR
C	    Define fputcfunc if fputc is a function rather than a subroutine
#	    ifdef hpux
#	        define fputcfunc
#	    endif
#	    ifdef absoft
#	        define fputcfunc
C	        Ensure that the GETFD function is defined for Absoft
	        INTEGER GETFD
#	    endif
#	    ifdef fputcfunc
C		On some systems, FPUTC is a function which
C		must be declared.
		INTEGER FPUTC
#	    endif
#	endif
	INTEGER IUNIT,NCH,NSP2,I
	CHARACTER CC,CH,STR*(*),OUTPUT*256
C
#	ifdef absoft
C	    The Absoft compiler does not generate code that automatically
C	    opens units until they are used in a Fortran I/O statement.
C	    Because we may be sending data for the first time using
C	    FPUTC (not normal Fortran I/O), we must explicitly open the
C	    unit if it isn't already opened.
	    IF (GETFD(IUNIT).LT.0) OPEN (IUNIT)
#	endif
C
C	Substitute spaces for "'s
	NCH = NSP2(STR)
	DO 10 I = 1,NCH
	    CH = STR(I:I)
	    IF (CH.EQ.'"') CH = ' '
	    OUTPUT(I:I) = CH
10	CONTINUE
C
#	ifdef unix
C	    On UNIX, carriage control may or may not work properly.  So
C	    it is emulated here using the system routine FPUTC.
	    IF (CC.EQ.'1') THEN
#		ifdef fputcfunc
C		    Actually FPUTC isn't all that standard.  On some
C		    systems (HP-UX) it is a function.
		    IERR = FPUTC(IUNIT,CHAR(12))
#		else
C		    On other systems (Linux) it is a subroutine.
		    CALL FPUTC(IUNIT,CHAR(12),IERR)
#		endif
	    ELSEIF (CC.EQ.'0') THEN
#		ifdef fputcfunc
		    IERR = FPUTC(IUNIT,CHAR(10))
#		else
		    CALL FPUTC(IUNIT,CHAR(10),IERR)
#		endif
	    ELSEIF (CC.EQ.'$') THEN
#		ifdef fputcfunc
		    IERR = FPUTC(IUNIT,CHAR(13))
#		else
		    CALL FPUTC(IUNIT,CHAR(13),IERR)
#		endif
	    ENDIF
	    DO 20 I = 1,NCH
#		ifdef fputcfunc
		    IERR = FPUTC(IUNIT,OUTPUT(I:I))
#		else
		    CALL FPUTC(IUNIT,OUTPUT(I:I),IERR)
#		endif
20	    CONTINUE
	    IF (CC.EQ.'+') THEN
#		ifdef fputcfunc
		    IERR = FPUTC(IUNIT,CHAR(13))
#		else
		    CALL FPUTC(IUNIT,CHAR(13),IERR)
#		endif
	    ELSEIF (CC.NE.'$'.AND.CC.NE.CHAR(0)) THEN
#		ifdef fputcfunc
		    IERR = FPUTC(IUNIT,CHAR(10))
#		else
		    CALL FPUTC(IUNIT,CHAR(10),IERR)
#		endif
	    ENDIF
#	else
C	    If not Unix, then we assume VMS where the output can be
C	    done with a simple write.
	    IF (NCH.LE.0) THEN
	        WRITE (IUNIT,20) CC
20	        FORMAT (A)
	    ELSE
	        WRITE (IUNIT,30) CC,OUTPUT(:NCH)
30	        FORMAT (A,A)
	    ENDIF
#	endif
	RETURN
	END
