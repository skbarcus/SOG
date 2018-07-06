      SUBROUTINE CPUBGN
C     Subroutine to return the elapsed CPU and clock times between
C     two points in a program in seconds.  To initialize, call CPUBGN.
C     At the end of the code being timed, call at the entry CPUEND
C     to obtain the CPU and clock time in that order.  This version
C     allows up to NCALL nested calls of CPUBGN.
C
C     Note that quite different methods of obtaining times are used
C     on VMS and Unix. Under Unix, it would be nice to find a way of
C     getting the clock time to a precision better than 1 second.
C
#     ifdef unix
C         Number of seconds per clock tick for ETIME routine (system
C         dependent)
#         ifdef hpux
              PARAMETER (CLKTIC = 1.)
#         else
              PARAMETER (CLKTIC = 1./60)
#         endif
#     endif
C     Number of nested calls allowed
      PARAMETER (NCALL = 10)
#     ifdef unix
C         Declarations specific to any Unix
          REAL SYS(NCALL),TARRAY(2)
          INTEGER TIME
#     else
C         Declarations specific to VMS
          CHARACTER TIMSTR*11
          INTEGER ITIME(2)
#     endif
      INTEGER IBUF(NCALL)
      DATA IPTR/0/
C
      IF (IPTR.GE.NCALL) RETURN
      IPTR = IPTR+1
#     ifdef unix
C         Get ellapsed cpu time
          SYS(IPTR) = ETIME(TARRAY)
C         Get system time in seconds
          IBUF(IPTR) = TIME()
#     else
C         On VMS use a system service
          IBUF(IPTR) = 0
          IRET = LIB$INIT_TIMER(IBUF(IPTR))
#     endif
      RETURN
C
      ENTRY CPUEND(TIMCPU,TIMCLK)
      IF (IPTR.LE.0) THEN
          TIMCPU = -1.
          TIMCLK = -1.
      ELSEIF (IPTR.GT.NCALL) THEN
          TIMCPU = -1.
          TIMCLK = -1.
          IPTR = IPTR-1
      ELSE
#         ifdef unix
              TIMCPU = (ETIME(TARRAY)-SYS(IPTR))*CLKTIC
              TIMCLK = TIME()-IBUF(IPTR)
#           else
C             If not Unix, must be VMS.
              IRET = LIB$STAT_TIMER(2,ITCPU,IBUF(IPTR))
              IRET = LIB$STAT_TIMER(1,ITIME,IBUF(IPTR))
              TIMCPU = ITCPU/100.
              IRET = SYS$ASCTIM(ITIML,TIMSTR,ITIME,%VAL(1))
              READ (TIMSTR,1) IHOUR,IMIN,SEC
1             FORMAT (I2,1X,I2,1X,F5.0)
              TIMCLK = (IHOUR*60.+IMIN)*60.+SEC
#         endif
          IPTR = IPTR-1
      ENDIF
      RETURN
      END
