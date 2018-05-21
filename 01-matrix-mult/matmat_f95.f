C **********************************************************************
C * Jonathan Senning <jonathan.senning@gordon.edu>                     *
C * Department of Mathematics and Computer Science                     *
C * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899           *
C *                                                                    *
C * FORTRAN 95 VERSION                                                 *
C *                                                                    *
C * Benchmark ijk, jki, and ikj matrix-matrix products.                *
C **********************************************************************
C
C $Smake: gfortran -std=f95 -O3 -o %F %f
C
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER N
      PARAMETER (N=500)
      DOUBLE PRECISION A(N,N),B(N,N)
      DOUBLE PRECISION C1(N,N),C2(N,N),C3(N,N)
      INTEGER ICMPPR
      INTEGER I,J
      REAL T1,T2
      REAL TIMIJK,TIMJKI,TIMIKJ
C
      WRITE (*,"(1X,'MATRIX-MATRIX MULTIPLY: MATRICES ARE ',
     +  I5,' X',I5)") N,N
C
C     INITIALIZE MATRICES
C
      DO J=1,N
         DO I=1,N
            CALL RANDOM_NUMBER(A(I,J))
            CALL RANDOM_NUMBER(B(I,J))
         END DO
      END DO
C
C     IJK PRODUCT
C
      CALL CPU_TIME(T1)
      CALL MATIJK(C1,A,B,N)
      CALL CPU_TIME(T2)
      TIMIJK=T2-T1
C
C     JKI PRODUCT
C
      CALL CPU_TIME(T1)
      CALL MATJKI(C2,A,B,N)
      CALL CPU_TIME(T2)
      TIMJKI=T2-T1
C
C     IKJ PRODUCT
C
      CALL CPU_TIME(T1)
      CALL MATIKJ(C3,A,B,N)
      CALL CPU_TIME(T2)
      TIMIKJ=T2-T1
C
C     OUTPUT RESULTS
C
      WRITE (*,*) '      IJK             JKI             IKJ'
      WRITE (*,*) '--------------  --------------  --------------'
      WRITE (*,"(1X,F10.6,' SEC',F12.6,' SEC',F12.6,' SEC')")
     +  TIMIJK,TIMJKI,TIMIKJ
      WRITE (*,"(1X,F10.2,' MFLOPS',F9.2,' MFLOPS',F9.2,' MFLOPS')")
     +  (2.0*N**3)/TIMIJK/1.0E6,(2.0*N**3)/TIMJKI/1.0E6,
     +  (2.0*N**3)/TIMIKJ/1.0E6
C
C     COMPARE PRODUCTS
C
      IF (ICMPPR(C1,C2,N).GT.0) THEN
         WRITE(*,"(1X,'C',I1,' /= C',I1,': VALIDATION ERROR')") 1,2
      ENDIF
      IF (ICMPPR(C1,C3,N).GT.0) THEN
         WRITE(*,"(1X,'C',I1,' /= C',I1,': VALIDATION ERROR')") 1,3
      ENDIF
      IF (ICMPPR(C2,C3,N).GT.0) THEN
         WRITE(*,"(1X,'C',I1,' /= C',I1,': VALIDATION ERROR')") 2,3
      ENDIF
C
C     ALL DONE
C
      END
C
C-----------------------------------------------------------------------------
C IJK MATRIX-MATRIX PRODUCT
C-----------------------------------------------------------------------------
C
      SUBROUTINE MATIJK(C,A,B,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),A(N,N),B(N,N)
      DO I=1,N
         DO J=1,N
            C(I,J)=0.0
            DO K=1,N
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
      RETURN
      END
C
C-----------------------------------------------------------------------------
C JKI MATRIX-MATRIX PRODUCT
C-----------------------------------------------------------------------------
C
      SUBROUTINE MATJKI(C,A,B,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),A(N,N),B(N,N)
      DO J=1,N
         DO I=1,N
            C(I,J)=0.0
         END DO
         DO K=1,N
            DO I=1,N   
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
      RETURN
      END
C
C-----------------------------------------------------------------------------
C IKJ MATRIX-MATRIX PRODUCT
C-----------------------------------------------------------------------------
C
      SUBROUTINE MATIKJ(C,A,B,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),A(N,N),B(N,N)
      DO I=1,N
         DO J=1,N
            C(I,J)=0.0
         END DO
         DO K=1,N
            DO J=1,N   
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO
      RETURN
      END
C
C-----------------------------------------------------------------------------
C COMPARE PRODUCTS
C-----------------------------------------------------------------------------
C
      INTEGER FUNCTION ICMPPR(C,D,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),D(N,N)
      DOUBLE PRECISION EPS
      INTEGER ICOUNT
      ICOUNT=0
      EPS=1.0E-12
      DO I=1,N
         DO J=1,N
            IF (ABS(C(I,J)-D(I,J))>EPS) THEN
               ICOUNT=ICOUNT+1
               WRITE (*,*) '****WARNING****'
               WRITE (*,"(1X,'(',I3,',',I3,'): ',F15.5,' <> ',F15.5)")
     +              I,J,C(I,J),D(I,J)
            ENDIF
         END DO
      END DO
      ICMPPR=ICOUNT
      END
