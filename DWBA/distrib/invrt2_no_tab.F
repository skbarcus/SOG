      SUBROUTINE INVRT2 (A,N,NDIM)
C     Subroutine to invert up to a NPMX X NPMX matrix. The matrix
C     is first scaled to avoid over(under)flow. A is the matrix
C     to invert.  When the routine returns, the inverse of A is
C     stored in A. N is the size of the matrix.  NDIM is the dimension
C     of A (A should be dimensioned A(NDIM,NDIM).  N may be less than
C     or equal to NDIM). This is the double precision version.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     The maximum matrix dimension that can be inverted
      PARAMETER (NPMX=100)
C
      DOUBLE PRECISION A(NDIM,NDIM),PIVOT(NPMX),DIAG(NPMX)
      INTEGER IPVOT(NPMX),INDEX(NPMX,2)
      DO 970 I=1,N
          DIAG(I)=ABS(A(I,I))
          IF (DIAG(I).EQ.0.D0) DIAG(I)=1.D0
          IPVOT(I)=0
970   CONTINUE
C     Scale the matrix using the diagonal elements
      DO 971 I=1,N
          DO 972 J=1,N
              A(J,I)=A(J,I)/SQRT(DIAG(I)*DIAG(J))
972       CONTINUE
971   CONTINUE
      DO 973 I=1,N
C         Search for pivot element
          T=0.D0
          DO 974 J=1,N
              IF (IPVOT(J).NE.1) THEN
                  DO 975 K=1,N
                      IF (IPVOT(K).LT.1.AND.ABS(T).LT.ABS(A(J,K))) THEN
                          IROW=J
                          ICOL=K
                          T=A(J,K)
                      ELSEIF (IPVOT(K).GT.1) THEN
                          GO TO 1
                      ENDIF
975               CONTINUE
              ENDIF
974       CONTINUE
          IPVOT(ICOL)=IPVOT(ICOL)+1
C         Put pivot element on diagonal
          IF (IROW.NE.ICOL) THEN
              DO 976 L=1,N
                  T=A(IROW,L)
                  A(IROW,L)=A(ICOL,L)
                  A(ICOL,L)=T
976           CONTINUE
          ENDIF
          INDEX(I,1)=IROW
          INDEX(I,2)=ICOL
          PIVOT(I)=A(ICOL,ICOL)
C         Divide pivot row by pivot element
          A(ICOL,ICOL)=1.D0
          DO 977 L=1,N
              A(ICOL,L)=A(ICOL,L)/PIVOT(I)
977       CONTINUE
C         Reduce non-pivot rows
          DO 978 M=1,N
              IF (M.NE.ICOL) THEN
                  T=A(M,ICOL)
                  A(M,ICOL)=0.D0
                  DO 979 L=1,N
                      A(M,L)=A(M,L)-A(ICOL,L)*T
979               CONTINUE
              ENDIF
978       CONTINUE
973   CONTINUE
C     Interchange columns
      DO 980 I=1,N
          L=N-I+1
          IF (INDEX(L,1).NE.INDEX(L,2)) THEN
              JROW=INDEX(L,1)
              JCOL=INDEX(L,2)
              DO 981 K=1,N
                  T=A(K,JROW)
                  A(K,JROW)=A(K,JCOL)
                  A(K,JCOL)=T
981           CONTINUE
          ENDIF
980   CONTINUE
C     Unscale the result using the diagonal elements
1     DO 982 I=1,N
          DO 983 J=1,N
              A(J,I)=A(J,I)/SQRT(DIAG(I)*DIAG(J))
983       CONTINUE
982   CONTINUE
      RETURN
      END
