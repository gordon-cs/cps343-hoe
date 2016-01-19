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
      DOUBLE PRECISION CS1,CS2,CS3
      DOUBLE PRECISION CMPPRD
      INTEGER I,J
      REAL T1,T2
      REAL TIMIJK,TIMIKJ,TIMJKI
      REAL TARRAY(2),ETIME
c      REAL RAND
C
      WRITE(*,5) N,N
   5  FORMAT(1X,'MATRIX-MATRIX MULTIPLY: MATRICES ARE ',I5,' X',I5)
C
C     INITIALIZE MATRICES
C
      DO 20 J=1,N
         DO 10 I=1,N
C            A(I,J)=0.1*MOD(2*I+5*J,10)
C            B(I,J)=0.1*MOD(4*I+3*J,10)
            A(I,J)=RAND(0)
            B(i,J)=RAND(0)
 10      CONTINUE
 20   CONTINUE
C
C     IJK PRODUCT
C
      T1 = ETIME(TARRAY)
      CALL MATIJK(C1,A,B,N)
      T2 = ETIME(TARRAY)
      TIMIJK=T2-T1
C
C     IKJ PRODUCT
C
      T1 = ETIME(TARRAY)
      CALL MATIKJ(C2,A,B,N)
      T2 = ETIME(TARRAY)
      TIMIKJ=T2-T1
C
C     JKI PRODUCT
C
      T1=ETIME(TARRAY)
      CALL MATJKI(C3,A,B,N)
      T2=ETIME(TARRAY)
      TIMJKI=T2-T1
C
C     OUTPUT RESULTS
C
      WRITE(*,30)
      WRITE(*,40)
      WRITE(*,50) TIMIJK,TIMJKI,TIMIKJ
      WRITE(*,60) (2.0*N**3)/TIMIJK/1.0E6,(2.0*N**3)/TIMJKI/1.0E6,
     *    (2.0*N**3)/TIMIKJ/1.0E6
  30  FORMAT(1X,'      IJK             JKI             IKJ')
  40  FORMAT(1X,'--------------  --------------  --------------')
  50  FORMAT(1X,F10.6,' SEC',F12.6,' SEC',F12.6, ' SEC')
  60  FORMAT(1X,F10.2,' MFLOPS',F9.2,' MFLOPS',F9.2,' MFLOPS')
C
C     COMPARE PRODUCTS
C
      CS1=CMPPRD(C1,C2,N)
      CS2=CMPPRD(C1,C3,N)
      CS3=CMPPRD(C2,C3,N)
      IF ((CS1.NE.CS2).OR.(CS1.NE.CS3).OR.(CS2.NE.CS3)) THEN
        WRITE(*,70) CS1,CS2,CS3
  70    FORMAT(1X,'CHECKSUM ERROR:',1X,F18.9,1X,F18.9,1X,F18.9)
      ENDIF
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
C COMPARE PRODUCTS
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION CMPPRD(C,D,N)
      INTEGER N
      DOUBLE PRECISION C(N,N),D(N,N)
      DOUBLE PRECISION CHECK
      CHECK=0.0
      DO 20 I=1,N
         DO 10 J=1,N
            CHECK=CHECK+C(I,J)
            IF (C(I,J).NE.D(I,J)) THEN
               WRITE(*,30)
               WRITE(*,40) I,J,C(I,J),D(I,J)
            ENDIF
 10      CONTINUE
 20   CONTINUE
      CMPPRD=CHECK
 30   FORMAT(1X,'****WARNING****')
 40   FORMAT(1X,'(',I3,',',I3,'): ',F15.5,' <> ',F15.5)
      RETURN
      END
