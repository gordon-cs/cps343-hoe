C **********************************************************************
C * Jonathan Senning <jonathan.senning@gordon.edu>                     *
C * Department of Mathematics and Computer Science                     *
C * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899           *
C *                                                                    *
C * FORTRAN 77 VERSION                                                 *
C *                                                                    *
C * Benchmark ijk, jki, and ikj matrix-matrix products.                *
C **********************************************************************
C
C $Smake: gfortran -std=legacy -O3 -o %F %f
C
C23456789012345678901234567890123456789012345678901234567890123456789012
      PROGRAM MAIN
      IMPLICIT CHARACTER*1 (A-Z)
      INTEGER N
      PARAMETER (N=500)
      DOUBLE PRECISION A(N,N),B(N,N)
      DOUBLE PRECISION C1(N,N),C2(N,N),C3(N,N)
      INTEGER ICMPPR
      INTEGER I,J
      INTEGER M,P,X
      REAL T1,T2,GFLPCT
      REAL TIMIJK,TIMJKI,TIMIKJ
      REAL TARRAY(2),ETIME
C
C     INITIALIZE MATRICES
C
      M=16807
      P=2147483647
      X=42
      DO 30 J=1,N
         DO 20 I=1,N
C            A(I,J)=RAND(0)
C            B(i,J)=RAND(0)
            X=MOD(M*X,P)
            A(I,J)=1.*X/P
            X=MOD(M*X,P)
            B(I,J)=1.*X/P
C            WRITE(*,32) A(I,J),B(I,J)
 20      CONTINUE
 30   CONTINUE
C 32   FORMAT(1X,'A=',F10.2,' B=',F10.2)

C
C     IJK PRODUCT
C
      T1 = ETIME(TARRAY)
      CALL MATIJK(C1,A,B,N)
      T2 = ETIME(TARRAY)
      TIMIJK=T2-T1
C
C     JKI PRODUCT
C
      T1=ETIME(TARRAY)
      CALL MATJKI(C2,A,B,N)
      T2=ETIME(TARRAY)
      TIMJKI=T2-T1
C
C     IKJ PRODUCT
C
      T1 = ETIME(TARRAY)
      CALL MATIKJ(C3,A,B,N)
      T2 = ETIME(TARRAY)
      TIMIKJ=T2-T1
C
C     OUTPUT RESULTS
C
      GFLPCT = 2.0*N**3/1.0E9
      WRITE(*,40) N,GFLPCT/TIMIJK,GFLPCT/TIMJKI,GFLPCT/TIMIKJ
 40   FORMAT(1X,'Fortran 77     (',i3,') ijk: ',F6.3,' gflops, jki: ',
     *  F6.3,' gflops, ikj: ',F6.3,' gflops')
C
C     COMPARE PRODUCTS
C
      IF (ICMPPR(C1,C2,N).GT.0) THEN
         WRITE(*,50) 1,2
      ENDIF
      IF (ICMPPR(C1,C3,N).GT.0) THEN
         WRITE(*,50) 1,3
      ENDIF
      IF (ICMPPR(C2,C3,N).GT.0) THEN
         WRITE(*,50) 2,3
      ENDIF
 50   FORMAT(1X,'C',I1,' /= C',I1,': VALIDATION ERROR')
C
C     ALL DONE
C
      STOP
      END
C
C-----------------------------------------------------------------------------
C IJK MATRIX-MATRIX PRODUCT
C-----------------------------------------------------------------------------
C
      SUBROUTINE MATIJK(C,A,B,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),A(N,N),B(N,N)
      DO 30 I=1,N
         DO 20 J=1,N
            C(I,J)=0.0
            DO 10 K=1,N
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
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
      DO 40 J=1,N
         DO 10 I=1,N
            C(I,J)=0.0
 10      CONTINUE
         DO 30 K=1,N
            DO 20 I=1,N
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
 20         CONTINUE
 30      CONTINUE
 40   CONTINUE
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
      DO 40 I=1,N
         DO 10 J=1,N
            C(I,J)=0.0
 10      CONTINUE
         DO 30 K=1,N
            DO 20 J=1,N
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
 20         CONTINUE
 30      CONTINUE
 40   CONTINUE
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
      DO 20 I=1,N
         DO 10 J=1,N
            IF (ABS(C(I,J)-D(I,J)).GT.EPS) THEN
               ICOUNT=ICOUNT+1
               WRITE(*,30)
               WRITE(*,40) I,J,C(I,J),D(I,J),C(I,J)-D(I,J)
            ENDIF
 10      CONTINUE
 20   CONTINUE
      ICMPPR=ICOUNT
 30   FORMAT(1X,'****ERROR****')
 40   FORMAT(1X,'(',I3,',',I3,'): ',F15.5,' /= ',F15.5,' Diff = ',E15.5)
      RETURN
      END
