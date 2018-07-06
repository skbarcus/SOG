C	This program converts the non-standard tabbed format of fortran
C	programs into the standard format in which the continuation
C	character is in column 6 and comands start in column 7.
c
C	Version 1B was slightly modified to make the code run on both
C	Unix and VMS.  Where differences are unavoidable, in-line comments
C	are used to indicate what should be compiled and what should be
C	commented out for particular systems.  At the moment the comments
C	can be:
C	    ! VMS       Only VMS
C	    ! Unix      Any Unix system including HP-UX (ie not VMS)
C	    ! NHPUX     Any Unix except HP-UX
C	The carriage control in format lines has been modified to produce
C	a consistent behavior on UNIX.  The standard "hello message"
C	is also printed.
C
C	Version 1C has been modified to use cpp conventions for conditional
C	compilation.  New standard subroutines for prompting and hello
C	message are installed.
C
	IMPLICIT INTEGER (A-Z)
C	Version number of this program
	CHARACTER VER*(*)
	PARAMETER (VER = '1C')
C	Last edit date for this program
	CHARACTER EDIT*(*)
	PARAMETER (EDIT = '10/10/96')
C	Output line lingth
	PARAMETER (OLEN = 80)
	CHARACTER INPUT*(OLEN),OUTPUT*(OLEN)
	LOGICAL TAB,CONTIN,STNCOM
C
C	Print a "Hello" message"
	CALL HELLOM(6,'NoTab',VER,EDIT,INPUT)
C
	CALL PROMPT(6,'Input file?')
	READ (5,11) INPUT
11	FORMAT (A)
#	ifdef linux
	    OPEN (UNIT=1,FILE=INPUT,STATUS='OLD')
#	else
	    OPEN (UNIT=1,FILE=INPUT,STATUS='OLD',READONLY)
#	endif
	CALL PROMPT(6,'Output file?')
	READ (5,11) INPUT
#	ifdef unix
	    OPEN (UNIT=2,FILE=INPUT,STATUS='UNKNOWN')
#	else
C	    We assume VMS if not Unix.
	    OPEN (UNIT=2,FILE=INPUT,STATUS='NEW',CARRIAGECONTROL='LIST',
	1    RECL=OLEN)
#	endif
C
C	Read next string form input file
1	READ (1,11,END=9999) INPUT
	IPTR = 0
	OPTR = 0
	ILEN = NSP2(INPUT)
	TAB = .FALSE.
	STNCOM = .FALSE.
	IF (ICHAR(INPUT(1:1)).EQ.12) THEN
C	    The line starts with a form feed.  Remove it.  If the line
C	    is then blank, forget the whole line.
	    IPTR = 1
	    IF (INPUT(2:ILEN).EQ.' ') GO TO 1
	ENDIF
30	IPTR = IPTR+1
	IF (IPTR.GT.ILEN) GO TO 100
	TAB = INPUT(IPTR:IPTR).EQ.CHAR(9)
	IF (INPUT(IPTR:IPTR).EQ.' '.OR.TAB) GO TO 200
	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = INPUT(IPTR:IPTR)
	STNCOM = .TRUE.
	GO TO 30
C
C	A space was seen
200	IPTR = IPTR+1
	IF (IPTR.GT.ILEN) GO TO 100
	IF (.NOT.TAB.AND.(INPUT(IPTR:IPTR).EQ.' '.AND.IPTR.LT.7))
	1 GO TO 200
C
C	A non-seperator after the first seperator was seen
210	CONTIN = (IPTR.EQ.6.AND..NOT.TAB).OR.(INPUT(IPTR:IPTR).GE.'0'
	1 .AND.INPUT(IPTR:IPTR).LE.'9'.AND.TAB).AND..NOT.STNCOM
	IF (CONTIN) GO TO 300
	COL = 0
410	IF (OPTR.GE.6) GO TO 400
	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = ' '
	GO TO 410
400	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = INPUT(IPTR:IPTR)
	IF (OUTPUT(OPTR:OPTR).NE.' '.AND.COL.EQ.0.AND..NOT.CONTIN)
	1 COL = OPTR
	IPTR = IPTR+1
	IF (IPTR.LE.ILEN) GO TO 400
C
C	Write out the output string
100	WRITE (2,11) OUTPUT(1:OPTR)
	GO TO 1
C
C	This is a continuation line.  Put it in the new format.
300	IF (OPTR.GE.5) GO TO 330
	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = ' '
	GO TO 300
330	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = '+'
340	IF (OPTR.GE.COL) GO TO 350
	OPTR = OPTR+1
	IF (OPTR.GT.OLEN) GO TO 100
	OUTPUT(OPTR:OPTR) = ' '
	GO TO 340
350	IPTR = IPTR+1
	IF (IPTR.GT.ILEN) GO TO 100
	IF (INPUT(IPTR:IPTR).EQ.' ') GO TO 350
	GO TO 400
C
9999	CLOSE (UNIT=1)
	CLOSE (UNIT=2)
	STOP 'All done'
	END
