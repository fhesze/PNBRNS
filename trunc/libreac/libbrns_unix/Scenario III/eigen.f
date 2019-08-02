c    $Id: eigen.f 16 2007-10-18 12:36:47Z centler $
c The following routines (based on the original EISPACK library)
c perform a diagonalization of a real symmetric matrix based on the
c QL algorithm.
cThe first step is to reduce the matrix to tridiagonal form (TRED2) and
cthen a further routine (TQLI) finds the eigenvalues and eigenvectors of
c     a tridiagonal matrix:
c     http://rsc.anu.edu.au/~harry/COURSES/MATHMETH/node70.html

c      PROGRAM D11R4
       subroutine eigen(c,np)
C     Driver for routine TQLI
c      PARAMETER(NP=10,TINY=1.0E-6)
       implicit none
       real*8 c,a,d,e,f,tiny
       integer i,j,np,k
       PARAMETER(TINY=3.d-6)
       DIMENSION A(NP,NP),C(NP,NP),D(NP),E(NP),F(NP)
c      DATA C/5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,-4.0,
c     *            4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0,
c     *            3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,-2.0,
c     *            2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,-1.0,
c     *            1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,0.0,
c     *            0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0,
c     *            -1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,
c     *            -2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,3.0,
c     *            -3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,4.0,
c     *            -4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0/
       DO 12 I=1,NP
          DO 11 J=1,NP
c            c(i,j)=0.
c            if(i.eq.j)c(i,j)=1.
             A(I,J)=C(I,J)
 11      CONTINUE
 12   CONTINUE

       CALL TRED2(A,NP,NP,D,E)
       CALL TQLI(D,E,NP,NP,A)
       WRITE(*,'(/1X,A)') 'Eigenvectors for a real symmetric matrix'
       DO 16 I=1,NP
          DO 14 J=1,NP
             F(J)=0.0d0
             DO 13 K=1,NP
                F(J)=F(J)+C(J,K)*A(K,I)
 13         CONTINUE
 14      CONTINUE
          WRITE(*,*) 'Eigenvalue',I,' =',D(I)
!         WRITE(*,'(/1X,A,I3,A,F10.6)') 'Eigenvalue',I,' =',D(I)
         WRITE(*,'(/1X,T7,A,T17,A,T31,A)') 'Vector','Mtrx*Vect.','Ratio'
          DO 15 J=1,NP
             IF (dabs(A(J,I)).LT.TINY) THEN
!            IF (ABS(A(J,I)).LT.TINY) THEN
                WRITE(*,'(1X,2F12.6,A12)') A(J,I),F(J),'div. by 0'
             ELSE
                WRITE(*,'(1X,2F12.6,E14.6)') A(J,I),F(J),
     *              F(J)/A(J,I)
             ENDIF
 15      CONTINUE
          WRITE(*,'(/1X,A)') 'press ENTER to continue...'
          READ(*,*)
 16   CONTINUE
       END


c     ==================================
       SUBROUTINE TRED2(A,N,NP,D,E)
       implicit none
       integer n,np,i,l,k,j
       real*8 a,d,e,h,scale,f,g,hh
       DIMENSION A(NP,NP),D(NP),E(NP)
       IF(N.GT.1)THEN
          DO 18 I=N,2,-1
             L=I-1
             H=0.d0
             SCALE=0.d0
             IF(L.GT.1)THEN
                DO 11 K=1,L
                   SCALE=SCALE+dabs(A(I,K))
!                  SCALE=SCALE+ABS(A(I,K))
 11            CONTINUE
                IF(SCALE.EQ.0.d0)THEN
                   E(I)=A(I,L)
                ELSE
                   DO 12 K=1,L
                      A(I,K)=A(I,K)/SCALE
                      H=H+A(I,K)**2
 12               CONTINUE
                   F=A(I,L)
                   G=-SIGN(dSQRT(H),F)
!                  G=-SIGN(SQRT(H),F)
                   E(I)=SCALE*G
                   H=H-F*G
                   A(I,L)=F-G
                   F=0.
                   DO 15 J=1,L
                      A(J,I)=A(I,J)/H
                      G=0.
                      DO 13 K=1,J
                         G=G+A(J,K)*A(I,K)
 13                  CONTINUE
                      IF(L.GT.J)THEN
                         DO 14 K=J+1,L
                            G=G+A(K,J)*A(I,K)
 14                     CONTINUE
                      ENDIF
                      E(J)=G/H
                      F=F+E(J)*A(I,J)
 15               CONTINUE
                   HH=F/(H+H)
                   DO 17 J=1,L
                      F=A(I,J)
                      G=E(J)-HH*F
                      E(J)=G
                      DO 16 K=1,J
                         A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
 16                  CONTINUE
 17               CONTINUE
                ENDIF
             ELSE
                E(I)=A(I,L)
             ENDIF
             D(I)=H
 18      CONTINUE
       ENDIF
       D(1)=0.d0
       E(1)=0.d0
       DO 23 I=1,N
          L=I-1
          IF(D(I).NE.0.d0)THEN
             DO 21 J=1,L
                G=0.
                DO 19 K=1,L
                   G=G+A(I,K)*A(K,J)
 19            CONTINUE
                DO 20 K=1,L
                   A(K,J)=A(K,J)-G*A(K,I)
 20            CONTINUE
 21         CONTINUE
          ENDIF
          D(I)=A(I,I)
          A(I,I)=1.d0
          IF(L.GE.1)THEN
             DO 22 J=1,L
                A(I,J)=0.d0
                A(J,I)=0.d0
 22         CONTINUE
          ENDIF
 23   CONTINUE
       RETURN
       END

       SUBROUTINE TQLI(D,E,N,NP,Z)
       implicit none
       integer itmax,i,n,np,m,l,iter,k
       real*8 d,e,z,dd,g,r,s,c,p,f,b
       DIMENSION D(NP),E(NP),Z(NP,NP)

c     max. number of iterations
       itmax=1000

       IF (N.GT.1) THEN
          DO 11 I=2,N
             E(I-1)=E(I)
 11      CONTINUE
          E(N)=0.
          DO 15 L=1,N
             ITER=0
 1          DO 12 M=L,N-1
                DD=dABS(D(M))+dABS(D(M+1))
                IF (dABS(E(M))+DD.EQ.DD) GO TO 2
!               DD=ABS(D(M))+ABS(D(M+1))
!               IF (ABS(E(M))+DD.EQ.DD) GO TO 2
 12         CONTINUE
             M=N
 2          IF(M.NE.L)THEN
c               IF(ITER.EQ.30)PAUSE 'too many iterations'
                IF(ITER.EQ.itmax)PAUSE 'too many iterations'
                ITER=ITER+1
                G=(D(L+1)-D(L))/(2.d0*E(L))
                R=dSQRT(G**2+1.)
!               R=SQRT(G**2+1.)
                G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
                S=1.d0
                C=1.d0
                P=0.d0
                DO 14 I=M-1,L,-1
                   F=S*E(I)
                   B=C*E(I)
                   IF(dABS(F).GE.dABS(G))THEN
!                  IF(ABS(F).GE.ABS(G))THEN
                      C=G/F
                      R=dSQRT(C**2+1.d0)
!                     R=SQRT(C**2+1.)
                      E(I+1)=F*R
                      S=1./R
                      C=C*S
                   ELSE
                      S=F/G
                      R=dSQRT(S**2+1.d0)
!                     R=SQRT(S**2+1.)
                      E(I+1)=G*R
                      C=1./R
                      S=S*C
                   ENDIF
                   G=D(I+1)-P
                   R=(D(I)-G)*S+2.d0*C*B
                   P=S*R
                   D(I+1)=G+P
                   G=C*R-B
                   DO 13 K=1,N
                      F=Z(K,I+1)
                      Z(K,I+1)=S*Z(K,I)+C*F
                      Z(K,I)=C*Z(K,I)-S*F
 13               CONTINUE
 14            CONTINUE
                D(L)=D(L)-P
                E(L)=G
                E(M)=0.d0
                GO TO 1
             ENDIF
 15      CONTINUE
       ENDIF
       RETURN
       END
