c    $Id: TRIDAG.F 16 2007-10-18 12:36:47Z centler $
c  ********************TRIDAG.F********************
       subroutine TRIDAG(A1,B1,C1,R1,U1,N)
 !     implicit double precision(a-h,o-z)
 !     PARAMETER (NMAX=10000)
 !     DIMENSION GAM(NMAX),A1(N),B1(N),C1(N),R1(N),U1(N)
 ! note the different array size declaration for the concentration

       implicit real*8 (a-h,o-z)

       DIMENSION GAM(N),A1(N),B1(N),C1(N),R1(N),U1(0:N+1)
!      DIMENSION GAM(N),A1(N),B1(N),C1(N),R1(N),U1(N)
       IF(B1(1).EQ.0.d0)PAUSE 'Crashed where b1(1)=0'
       BET=B1(1)
       U1(1)=R1(1)/BET
       DO 11 J=3,N,2
         GAM(J)=C1(J-2)/BET
         BET=B1(J)-A1(J)*GAM(J)
         IF(BET.EQ.0.d0)PAUSE 'Crashed where BET=0'
         U1(J)=(R1(J)-A1(J)*U1(J-2))/BET
11    CONTINUE
       DO 12 J=N-2,1,-2
         U1(J)=U1(J)-GAM(J+2)*U1(J+2)
12    CONTINUE
       RETURN
       END
c********************************************************
