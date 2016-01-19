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
      REAL T1,T2
      REAL TIMIJK,TIMJKI,TIMIKJ
      REAL TARRAY(2),ETIME
c      REAL RAND
C
      WRITE(*,10) N,N
 10   FORMAT(1X,'MATRIX-MATRIX MULTIPLY: MATRICES ARE ',I5,' X',I5)
C
C     INITIALIZE MATRICES
C
      DO 30 J=1,N
         DO 20 I=1,N
            A(I,J)=RAND(0)
            B(i,J)=RAND(0)
 20      CONTINUE
 30   CONTINUE
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
      WRITE(*,40)
      WRITE(*,50)
      WRITE(*,60) TIMIJK,TIMJKI,TIMIKJ
      WRITE(*,70) (2.0*N**3)/TIMIJK/1.0E6,(2.0*N**3)/TIMJKI/1.0E6,
     *     (2.0*N**3)/TIMIKJ/1.0E6
  40  FORMAT(1X,'      IJK             JKI             IKJ')
  50  FORMAT(1X,'--------------  --------------  --------------')
  60  FORMAT(1X,F10.6,' SEC',F12.6,' SEC',F12.6, ' SEC')
  70  FORMAT(1X,F10.2,' MFLOPS',F9.2,' MFLOPS',F9.2,' MFLOPS')
C
C     COMPARE PRODUCTS
C
      IF (ICMPPR(C1,C2,N).GT.0) THEN
         WRITE(*,80) 1,2
      ENDIF
      IF (ICMPPR(C1,C3,N).GT.0) THEN
         WRITE(*,80) 1,3
      ENDIF
      IF (ICMPPR(C2,C3,N).GT.0) THEN
         WRITE(*,80) 2,3
      ENDIF
 80   FORMAT(1X,'VALIDATION ERROR: C',I1,' /= C',I1)
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
      INTEGER ICOUNT
      ICOUNT=0
      DO 20 I=1,N
         DO 10 J=1,N
            IF (C(I,J).NE.D(I,J)) THEN
               ICOUNT=ICOUNT+1
               WRITE(*,30)
               WRITE(*,40) I,J,C(I,J),D(I,J)
            ENDIF
 10      CONTINUE
 20   CONTINUE
      ICMPPR=ICOUNT
 30   FORMAT(1X,'****ERROR****')
 40   FORMAT(1X,'(',I3,',',I3,'): ',F15.5,' /= ',F15.5)
      RETURN
      END
