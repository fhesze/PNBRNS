c    $Id: frngdb.f 16 2007-10-18 12:36:47Z centler $
       SUBROUTINE RNFARR(AA,N)
C       FORTRAN 77 version of "ranf_array"
C       from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
C       including the MODIFICATIONS made in the 9th printing (2002)
C       ********* see the book for explanations and caveats! *********
       IMPLICIT DOUBLE PRECISION (A,R,Y)
       DIMENSION AA(*)
       PARAMETER (KK=100)
       PARAMETER (LL=37)
       COMMON /RSTATE/ RANX(KK)
       SAVE /RSTATE/
       DO 1 J=1,KK
 1    AA(J)=RANX(J)
       DO 2 J=KK+1,N
          Y=AA(J-KK)+AA(J-LL)
          AA(J)=Y-IDINT(Y)
 2    CONTINUE
       DO 3 J=1,LL
          Y=AA(N+J-KK)+AA(N+J-LL)
          RANX(J)=Y-IDINT(Y)
 3    CONTINUE
       DO 4 J=LL+1,KK
          Y=AA(N+J-KK)+RANX(J-LL)
          RANX(J)=Y-IDINT(Y)
 4    CONTINUE
       END

       SUBROUTINE RNFSTR(SEED)
       IMPLICIT DOUBLE PRECISION (A,R,U,V)
       IMPLICIT INTEGER (T)
       PARAMETER (KK=100)
       PARAMETER (LL=37)
       PARAMETER (MM=2**30)
       PARAMETER (ULP=1D0/(2D0**52))
       PARAMETER (TT=70)
       PARAMETER (KKK=KK+KK-1)
       INTEGER SEED,S,SSEED
       DOUBLE PRECISION SS
       DIMENSION U(KKK)
       COMMON /RSTATE/ RANX(KK)
       SAVE /RSTATE/
       IF (SEED .LT. 0) THEN
          SSEED=MM-1-MOD(-1-SEED,MM)
       ELSE
          SSEED=MOD(SEED,MM)
       END IF
       SS=2D0*ULP*(SSEED+2)
       DO 1 J=1,KK
          U(J)=SS
          SS=SS+SS
          IF (SS .GE. 1D0) SS=SS-1D0+2*ULP
 1    CONTINUE
       U(2)=U(2)+ULP
       S=SSEED
       T=TT-1
 10   DO 12 J=KK,2,-1
          U(J+J-1)=U(J)
 12      U(J+J-2)=0
       DO 14 J=KKK,KK+1,-1
          V=U(J-(KK-LL))+U(J)
          U(J-(KK-LL))=V-IDINT(V)
          V=U(J-KK)+U(J)
          U(J-KK)=V-IDINT(V)
 14   CONTINUE
       IF (MOD(S,2) .EQ. 1) THEN
          DO 16 J=KK,1,-1
 16         U(J+1)=U(J)
          U(1)=U(KK+1)
          V=U(LL+1)+U(KK+1)
          U(LL+1)=V-IDINT(V)
       END IF
       IF (S .NE. 0) THEN
          S=S/2
       ELSE
          T=T-1
       END IF
       IF (T .GT. 0) GO TO 10
       DO 20 J=1,LL
 20      RANX(J+KK-LL)=U(J)
       DO 21 J=LL+1,KK
 21      RANX(J-LL)=U(J)
       DO 22 J=1,10
 22      CALL RNFARR(U,KKK)
       END

