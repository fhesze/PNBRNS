c    $Id: NEWT.F 16 2007-10-18 12:36:47Z centler $
       SUBROUTINE newt(check,jloc,iflag,tolf,tolx,newton,ilud)
      include 'common_geo.inc'
      include 'common.inc'

c***********************************************************************
c      NEWT: newton iteration with linesearch                          *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - common blocks added, so variables don't get undefined         *
c      - NMAX replaced by ncomp                                        *
c      - tolerance levels and iteration numbers from global definitons *
c      - variable declarations changed from REAL to REAL*8             *
c      - preserving fixed concentration boundaries (lnsrch)            *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************
      integer ilud
      integer jloc,newton,iflag
       LOGICAL check

       REAL*8 fvec,TOLF,TOLMIN,TOLX,STPMX
!x      PARAMETER (TOLMIN=1.0d-6,STPMX=100.0d0)
       PARAMETER (TOLMIN=1.d-12,STPMX=100.0d0)
       COMMON /newtv/ fvec(ncomp)
       SAVE /newtv/
CU    USES fdjac,fmin,lnsrch,lubksb,ludcmp
       INTEGER i,its,j,indx(ncomp)
       REAL*8 d,den,f,fold,stpmax,sum,temp,test,fjac(ncomp,ncomp)
       REAL*8 g(ncomp),p(ncomp),xold(ncomp),fmin
       EXTERNAL fmin
       f=fmin(jloc)

       test=0.0d0
       do 11 i=1,ncomp
         if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
11    continue
       if(test.lt.0.01d0*TOLF)then
         check=.false.
         return
       endif
       sum=0.0d0
       do 12 i=1,ncomp
         sum=sum+sp(i,jloc)*sp(i,jloc)
12    continue
       stpmax=STPMX*max(dsqrt(sum),dfloat(ncomp))
       do 21 its=1,newton
        call jacobian(fjac,jloc)
         do 14 i=1,ncomp
           sum=0.0d0
           do 13 j=1,ncomp
             sum=sum+fjac(j,i)*fvec(j)
13        continue
           g(i)=sum
14      continue
         do 15 i=1,ncomp
           xold(i)=sp(i,jloc)
15      continue
         fold=f
         do 16 i=1,ncomp
           p(i)=-fvec(i)
16      continue

c     LUD choice
      if(ilud.eq.1)then
           call ludcmp(fjac,ncomp,ncomp,indx,d)
           call lubksb(fjac,ncomp,ncomp,indx,p)
      else
           call DGESV( ncomp, 1, fjac, ncomp, indx, p, ncomp, INFO )
           if(info.ne.0)then
                if(info.lt.0)write(*,*)'illegal value in LUD @:', info
                if(info.gt.0)write(*,*)'singular matrix in LUD'
           endif
      endif

c       do i=1,ncomp
c         if (p(i).lt.0.0d0) then
c           write(*,*) 'conc < 0 => reset to old', p(i)
c           p(i) = xold(i)
c         end if
c        end do

         call lnsrch(xold,fold,g,p,f,stpmax,check,fmin,jloc,tolx)
         test=0.d0 
         do 17 i=1,ncomp

! check for neg. conc!!!!!!
          ! avoid looping due to fixed conc boundary
          if ((jloc.eq.1) .and.(ibc(i,1).eq.0)) fvec(i) = 0.0d0
          if ((jloc.eq.nx).and.(ibc(i,2).eq.0)) fvec(i) = 0.0d0

           if(dabs(fvec(i)).gt.test)test=dabs(fvec(i))
17      continue
         if(test.lt.TOLF)then
           check=.false.
           return
         endif
         if(check)then
           test=0.d0
           den=max(f,0.5d0*ncomp)
           do 18 i=1,ncomp
             temp=dabs(g(i))*max(dabs(sp(i,jloc)),1.0d0)/den
             if(temp.gt.test)test=temp
18        continue
           if(test.lt.TOLMIN)then
             check=.true.
           else
             check=.false.
           endif
           return
         endif
         test=0.0d0
         do 19 i=1,ncomp
           temp=(dabs(sp(i,jloc)-xold(i)))/max(dabs(sp(i,jloc)),1.0d0)

          ! avoid looping due to fixed conc boundary
          if ((jloc.eq.1) .and.(ibc(i,1).eq.0)) temp = 0.0d0
          if ((jloc.eq.nx).and.(ibc(i,2).eq.0)) temp = 0.0d0
          if (sp(i,jloc).lt.0.0d0) iflag = 1 

           if(temp.gt.test)test=temp
19      continue
         if (test.lt.TOLX) return
21    continue
       ! flag if exceeding newton iterations
      if (iflag.eq.0) iflag = 2
      if (iflag.eq.1) iflag = 3
       END


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       SUBROUTINE lnsrch(xold,fold,g,p,f,stpmax,check,func,jloc,tolx)

c***********************************************************************
c      LNSRCH: Part of NEWT                                            *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - common blocks added, so variables don't get undefined         *
c      - NMAX replaced by np, which is equal to ncomp                  *
c      - variable declarations changed from REAL to REAL*8             *
c      - preserving fixed concentration boundaries                     *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************

         include 'common_geo.inc'
         include 'common.inc'

         INTEGER jloc
         LOGICAL check
        REAL*8 f,fold,stpmax,g(ncomp),p(ncomp),xold(ncomp),func,ALF,TOLX
!x        PARAMETER (ALF=1.e-4)
         PARAMETER (ALF=1.d-4)
         EXTERNAL func
CU    USES func
         INTEGER i
         REAL*8 a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum
         REAL*8 temp,test,tmplam

         check=.false.
         sum=0.0d0
         do 11 i=1,ncomp
           sum=sum+p(i)*p(i)
 11     continue
         if (sum.ge.0.0d0) then
          sum=dsqrt(sum)
        else
          write(*,*) 'sum < 0', sum
        end if
         if(sum.gt.stpmax)then
           do 12 i=1,ncomp
             p(i)=p(i)*stpmax/sum
 12       continue
         endif
         slope=0.0d0
         do 13 i=1,ncomp
           slope=slope+g(i)*p(i)
 13     continue
         test=0.0d0
         do 14 i=1,ncomp
           temp=dabs(p(i))/max(dabs(xold(i)),1.0d0)

!         if fixed concentration preserve boundary condition
           if ((jloc.eq.1) .and.(ibc(i,1).eq.0)) temp=0.0d0
           if ((jloc.eq.nx).and.(ibc(i,2).eq.0)) temp=0.0d0

           if(temp.gt.test)test=temp
 14     continue
         if (test.gt.0.0d0) then
          alamin=TOLX/test
        else
          alamin = 10.0d-12
        end if
         alam=1.0d0
 1     continue
         do 15 i=1,ncomp
           sp(i,jloc)=xold(i)+alam*p(i)

!         if fixed concentration preserve boundary condition
           if ((jloc.eq.1) .and.(ibc(i,1).eq.0)) sp(i,jloc) = spb(i,1)
           if ((jloc.eq.nx).and.(ibc(i,2).eq.0)) sp(i,jloc) = spb(i,2)

15      continue 
         f=func(jloc)
         if(alam.lt.alamin)then
           do 16 i=1,ncomp
             sp(i,jloc)=xold(i)
!           if fixed concentration preserve boundary condition
             if ((jloc.eq.1) .and.(ibc(i,1).eq.0)) sp(i,jloc) = spb(i,1)
             if ((jloc.eq.nx).and.(ibc(i,2).eq.0)) sp(i,jloc) = spb(i,2)
16        continue
           check=.true.
           return
         else if(f.le.fold+ALF*alam*slope)then
           return
         else
           if(alam.eq.1.)then
             tmplam=-slope/(2.0d0*(f-fold-slope))
           else
             rhs1=f-fold-alam*slope
             rhs2=f2-fold2-alam2*slope
c            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
c            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
             a=(rhs1/alam/alam-rhs2/alam2/alam)/(alam-alam2)
             b=(-alam2*rhs1/alam/alam+alam*rhs2/alam2/alam)/(alam-alam2)
             if(a.eq.0.0d0)then
               tmplam=-slope/(2.0d0*b)
             else
               disc=b*b-3.0d0*a*slope
               if(disc.lt.0.0d0) then
                 write(*,*) 'roundoff problem in lnsrch',jloc
                disc = 0.0d0
               end if
               tmplam=(-b+dsqrt(disc))/(3.0d0*a)
             endif
             if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
           endif
         endif
         alam2=alam
         f2=f
         fold2=fold
         alam=max(tmplam,0.1d0*alam)
         goto 1

         END

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

       FUNCTION fmin(jloc)

c**********************************************************************
c    FMIN: part of newt (newton&linesearch)                           *
c    from Press et al., documented there                              *
c    changes made:                                                    *
c    - added common block to preserve variables                       *
c    - changed calling statement and content to common block "newtv"  *
c      due to known variable in common block                          *
c    - promoted REAL to REAL*8                                        *
c**********************************************************************

       include 'common_geo.inc'
       include 'common.inc'

       INTEGER jloc
       REAL*8 fmin,fvec
       COMMON /newtv/ fvec(ncomp)
       SAVE /newtv/
CU    USES funcv
       INTEGER i
       REAL*8 sum

       call residual(fvec,jloc)
       sum=0.0d0
       do 11 i=1,ncomp
         sum=sum+fvec(i)*fvec(i)
11    continue
       fmin=0.5d0*sum

       return
       END

! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
