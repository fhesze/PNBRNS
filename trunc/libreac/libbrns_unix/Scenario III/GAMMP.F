c    $Id: GAMMP.F 16 2007-10-18 12:36:47Z centler $
! collection of functions and subroutines from Press et al.
! to calculate probabilities related to chi^2
! changes: 
!  pause -> write
!  real -> real*8 (& functions -> dfunctions)
!  loop labels
! CM, Dec 2002

       FUNCTION gammp(a,x) ! a is half the degrees of
! freedom, x is half of chi^2
       REAL*8 a,gammp,x
CU    USES gcf,gser
       REAL*8 gammcf,gamser,gln
!      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'

       if(x.lt.0.d0.or.a.le.0.d0)write(*,*) 'bad arguments in gammp'
       if(x.lt.a+1.d0)then

         call gser(gamser,a,x,gln)
         gammp=gamser
       else
         call gcf(gammcf,a,x,gln)
         gammp=1.-gammcf
       endif
       return
       END


       SUBROUTINE gcf(gammcf,a,x,gln)
       INTEGER ITMAX
       REAL*8 a,gammcf,gln,x,EPS,FPMIN
       PARAMETER (ITMAX=100,EPS=3.d-16,FPMIN=1.d-30)
!x      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
       INTEGER i
       REAL*8 an,b,c,d,del,h,gammln
       gln=gammln(a)
       b=x+1.d0-a
       c=1.d0/FPMIN
       d=1.d0/b
       h=d
       do 11 i=1,ITMAX
         an=-i*(i-a)
         b=b+2.d0
         d=an*d+b
         if(dabs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(dabs(c).lt.FPMIN)c=FPMIN
         d=1./d
         del=d*c
         h=h*del
         if(dabs(del-1.).lt.EPS)goto 1
11    continue
!      pause 'a too large, ITMAX too small in gcf'
       write(*,*) 'a too large, ITMAX too small in gcf/gammp'
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
       return
       END


       SUBROUTINE gser(gamser,a,x,gln)
       INTEGER ITMAX
       REAL*8 a,gamser,gln,x,EPS
       PARAMETER (ITMAX=100,EPS=3.d-16)
!x      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
       INTEGER n
       REAL*8 ap,del,sum,gammln
       gln=gammln(a)
       if(x.le.0.d0)then
!        if(x.lt.0.)pause 'x < 0 in gser'
         if(x.lt.0.d0)write(*,*) 'x < 0 in gser/gammp'
         gamser=0.d0
         return
       endif
       ap=a
       sum=1.d0/a
       del=sum
       do 22 n=1,ITMAX
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(dabs(del).lt.dabs(sum)*EPS)goto 2
22    continue
!      pause 'a too large, ITMAX too small in gser'
       write(*,*) 'a too large, ITMAX too small in gser/gammp'
2     gamser=sum*dexp(-x+a*dlog(x)-gln)
       return
       END



       FUNCTION gammln(xx)
       REAL*8 gammln,xx
       INTEGER j
!      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
       real*8 ser,stp,tmp,x,y,cof(6)
       SAVE cof,stp
       DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
       x=xx
       y=x
       tmp=x+5.5d0
       tmp=(x+0.5d0)*dlog(tmp)-tmp
       ser=1.000000000190015d0
       do 33 j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
33    continue
       gammln=tmp+dlog(stp*ser/x)
       return
       END
