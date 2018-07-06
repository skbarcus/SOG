      SUBROUTINE TOD(TIM,DAY)
C     Subroutine to return the day and date in standard
C     (American) format in two 8-character strings
C
C     Modified to compile with gfortran
C
      IMPLICIT NONE
      CHARACTER TIM*(*),DAY*(*),IARAY*12,TIMI*30
#     ifdef gfortran
          INTEGER*8 ITIM
          INTEGER IDAT(3)
          EQUIVALENCE (DA,IDAT(1)),(MN,IDAT(2)),(YR,IDAT(3))
#     endif
      INTEGER MN,DA,YR
C
#     ifdef gfortran
          ITIM = TIME8()
          CALL CTIME(ITIM,TIMI)
          TIM = TIMI(12:19)
          CALL IDATE(IDAT)
#else
          CALL TIME(TIMI)
          TIM = TIMI
          CALL IDATE(MN,DA,YR)
#endif
      WRITE (IARAY,1) MN+100,DA+100,YR+100
1     FORMAT (3I4)
      DAY = ' '
      DAY(1:2) = IARAY(3:4)
      DAY(3:3) = '/'
      DAY(4:5) = IARAY(7:8)
      DAY(6:6) = '/'
      DAY(7:8) = IARAY(11:12)
      RETURN
      END
