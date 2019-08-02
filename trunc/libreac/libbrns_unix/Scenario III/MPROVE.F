c    $Id: MPROVE.F 16 2007-10-18 12:36:47Z centler $
       SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
!      implicit double precision (a-h,o-z)

!      PARAMETER (NMAX=100)
!      DIMENSION A(NP,NP),ALUD(NP,NP),INDX(N),B(N),X(N),R(NMAX)

       implicit real*8 (a-h,o-z)

       DIMENSION A(NP,NP),ALUD(NP,NP),INDX(N),B(N),X(N),R(N)
c     REAL*8 SDP
       DO 12 I=1,N
         SDP=-B(I)
         DO 11 J=1,N
           SDP=SDP+DBLE(A(I,J))*DBLE(X(J))
11      CONTINUE
         R(I)=SDP
12    CONTINUE
       CALL LUBKSB(ALUD,N,NP,INDX,R)
       DO 13 I=1,N
         X(I)=X(I)-R(I)
13    CONTINUE
       RETURN
       END
c*****************************************************************
