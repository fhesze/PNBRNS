c    $Id: optimsa.f 16 2007-10-18 12:36:47Z centler $
      subroutine optimsa()

c***********************************************************************
c      OPTIMSA: Drives optimization using simulated annealing (AMOESBA)*
c                                                                      *
c      CM, Jan 2002                                                    *
c***********************************************************************

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

      integer iter, i, j, k, nrmeastot,itot,ncool
      real*8 funk,p(nopt+1,nopt),psum(nopt),y(nopt+1),yb,pb(nopt)
      real*8 alphatemp,temptr

      external funk

c***********************************************************************
c      initialize the vertices of the simplex                          *
c      the initial profile is taken from maple input                   *
c      parameters to be optimized are selected and put into matrix p   *
c      to create the vertices, parameters are perturbed by a factor 2  *
c      should the parameter be set to 0, it is arbitrarily set to 1    *
c       ...could be improved, I guess                                  *
c      ntotparam: total number of parameters                           *
c      idpar: index of parameter that gets optimized                   *
c      par: vector containing all parameters                           *
c      p: matrix, containing the list of optimizable parameters (rows) *
c         the 1st row are original parameters                          *
c         the 2nd row has a perturbated 1st parameter                  *
c         the 3rd row has a perturbated 2nd parameter a.s.o.           *
c***********************************************************************
        do k=1,ntotparam
         do i=1,nopt
            if (idpar(i).eq.k) then
              do j=1,nopt+1
               p(j,i) = par(k)
                   if ((i.eq.j-1).and.(j.gt.1)) then
                 p(j,i) = p(j,i)*(1.d0+perturb)
                  if (p(j,i).eq.0.d0) p(j,i) = 1.d0
                end if
             end do
           end if
         end do
        end do

c***********************************************************************
c      calculated obj. function associated with the initial vertices   *
c***********************************************************************
       do j=1, nopt+1
         do i=1, nopt
           psum(i) = p(j,i)
         end do
         y(j) = funk(psum)
      end do

c***********************************************************************
c      optimization using simualted annealing                          *
c      define starting temperature, and optimize at this temperature   *
c      for a fraction of the total allowed iterations (->itmax*1/ncool)*
c      then decrease temp, thus discourage uphill moves and continue   *
c      the starting temperature is estimated to allow ~80% uphill moves*
c      -> sum_nopt(1/nopt * exp(-deltaOF/temp0)) = 0.8                 *
c         truncated taylor series expansion gives:                     *
c         temp0 = 5*av_deltaOF, where temp0 is the starting temperature*
c         http://petaxp.rug.ac.be/~erik/research/research-part3.html   *
c      yb: initialized high, coming out contains lowest obj. function  *
c      pb: contains coming out the best parameter set                  *
c      ncool: number of  cooling steps                                 *
c      alphatemp: cooling speed. the larger it is, the more time the   *
c                 algorithm spends at low temperatures                 *
c      temptr: (current) annealing temperature                         *
c      iter: in = number of function evaluations at one temperature    *
c            if iter > 0 when coming out then coinvergence achieved    *
c      ittot: iterations performed                                     *
c***********************************************************************
       yb = 999999999.d0
      ittot = 0
      alphatemp = 2.d0
      ncool = 10
      temptr = 0.d0
      do i=2,nopt+1
        temptr = temptr + dabs(y(i) - y(1))
      end do
       temptr = 5.d0*temptr/nopt
       DO WHILE (ittot.lt.itmax)
!      do i=1,ncool
        iter = nint(dfloat(itmax)/dfloat(ncool))
        ittot = ittot + iter
c     write(3,*) 'temptr', temptr,ittot,iter,itmax
      write(*,*) 'temptr', temptr,ittot,iter,itmax
         call amebsa(p,y,nopt+1,nopt,nopt,pb,yb,funk,iter,temptr)
c     write(3,*) 'temptr', temptr,ittot,iter,itmax
      write(*,*) 'temptr', temptr,ittot,iter,itmax
        temptr = temptr*(1.d0-dfloat(ittot)/dfloat(itmax))**alphatemp
        if(iter.ge.0) ittot = itmax
       end do

c***********************************************************************
c      put the optimized parameters back into the complete list        *
c      assign them to their real names (transferback)                  *
c      and give final output                                           *
c***********************************************************************
        do k=1,ntotparam
          do i=1,nopt
            if (idpar(i).eq.k) then
!            par(k) = p(1,i)
             par(k) = pb(i)
              end if
          end do
        end do
        call transferback()

      write (*,*) 'BEST PARAMETER VALUES, OF:', yb
      write (3,*) 'BEST PARAMETER VALUES, OF:', yb
      do i=1,nopt
        write (*,*) i, par(idpar(i))
        write (3,*) i, par(idpar(i))
      end do

      return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       SUBROUTINE amebsa(p,y,mp,np,ndim,pb,yb,funk,iter,temptr)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

c***********************************************************************
c      AMEBSA: Parameter optimization using a simulated annealing      *
c      This subroutine is from Press et al. and documented there       *
c      The following changes were made:                                *
c      - common blocks added, so variables don't get undefined         *
c      - NMAX replaced by np, which is equal to ncomp                  *
c      - variable declarations changed from REAL to REAL*8             *
c      - renamed ftol to optftol, to avoid conflict w. common block    *
c      - ITMAX from common block common_drive                          *
c      - if exceeding ITMAX return instead of stop                     *
c      - changing stopping criterion (one added based on chi^2)        *
c                                                                      *
c                                                                      *
c      CM, Jan 2001                                                    *
c***********************************************************************

       INTEGER iter,mp,ndim,np,t
       REAL*8 temptr,yb,p(mp,np),pb(np),y(mp),funk
       EXTERNAL funk
CU    USES amotsa,funk,ran1
       INTEGER i,idum,ihi,ilo,inhi,j,m,n
       REAL*8 rtol,sum,swap,tt,yhi,ylo,ynhi,ysave,yt,ytry,psum(np),
     *amotsa,ran1
       COMMON /ambsa/ tt,idum

c      some initial issues
        if (iof.eq.1) then
          idegf = 0
          do t=1,ntopt
           do k=1,nrspmeas(t)
              do i=1,nrxmeas(t,k)
               idegf=idegf+1
             end do
           end do
         end do
          idegf = idegf - nopt
       else
         prob = 2.d0 ! any value greater than 1, so it doesn't interfere
       end if

       tt=-temptr
1     do 12 n=1,ndim
         sum=0.d0
         do 11 m=1,ndim+1
           sum=sum+p(m,n)
11      continue
         psum(n)=sum
12    continue
2     ilo=1
       inhi=1
       ihi=2
       ylo=y(1)+tt*dlog(ran1(idum))
       ynhi=ylo
       yhi=y(2)+tt*dlog(ran1(idum))
       if (ylo.gt.yhi) then
         ihi=1
         inhi=2
         ilo=2
         ynhi=yhi
         yhi=ylo
         ylo=ynhi
       endif

       do 13 i=3,ndim+1
         yt=y(i)+tt*dlog(ran1(idum))
         if(yt.le.ylo) then
           ilo=i
           ylo=yt
         endif
         if(yt.gt.yhi) then
           inhi=ihi
           ynhi=yhi
           ihi=i
           yhi=yt
         else if(yt.gt.ynhi) then
           inhi=i
           ynhi=yt
         endif
13    continue

       rtol=2.d0*dabs(yhi-ylo)/(dabs(yhi)+dabs(ylo))
      if (iof.eq.1) prob = gammp(dfloat(idegf)/2.d0,y(ilo)/2.d0)
       if((rtol.lt.optftol).or.(prob.lt.problim).or.(iter.lt.0)) then
!      if (rtol.lt.ftol.or.iter.lt.0) then
         if (rtol.lt.optftol) write(3,*) 'rtol<optftol',rtol,optftol
         if (prob.lt.problim) write(3,*) 'prob<problim',prob,problim
         if (iter.lt.0) write(3,*) 'iter<0'
         swap=y(1)
         y(1)=y(ilo)
         y(ilo)=swap
         do 14 n=1,ndim
           swap=p(1,n)
           p(1,n)=p(ilo,n)
           p(ilo,n)=swap
          if(yb.gt.y(ilo)) pb(n) = p(ilo,n)
14      continue
        ! in case solution was found in inital vertices
         if (yb.gt.y(ilo)) yb = y(ilo)
         return
       endif

       iter=iter-2
       ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,-1.0d0)
       if (ytry.le.ylo) then
         ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,2.0d0)
       else if (ytry.ge.ynhi) then
         ysave=yhi
         ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,0.5d0)
         if (ytry.ge.ysave) then
           do 16 i=1,ndim+1
             if(i.ne.ilo)then
               do 15 j=1,ndim
                 psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                 p(i,j)=psum(j)
15            continue
               y(i)=funk(psum)
             endif
16        continue
           iter=iter-ndim
           goto 1
         endif
       else
         iter=iter+1
       endif
       goto 2

       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       FUNCTION amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,fac)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

c***********************************************************************
c      AMOTSA: Part of AMEBSA (simulated annealing)                    *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - common blocks added, so variables don't get undefined         *
c      - NMAX replaced by np, which is equal to ncomp                  *
c      - variable declarations changed from REAL to REAL*8             *
c                                                                      *
c      CM, Jan 2002                                                    *
c***********************************************************************

       INTEGER ihi,mp,ndim,np,idum,j
       REAL*8 amotsa,fac,yb,yhi,p(mp,np),pb(np),psum(np),y(mp),funk
       EXTERNAL funk
CU    USES funk,ran1
       REAL*8 fac1,fac2,tt,yflu,ytry,ptry(np),ran1
       COMMON /ambsa/ tt,idum
      real*8 facmin,facmax,small,fact                      ! addition
      integer iwarn                                        ! addition

!      fac1=(1.-fac)/ndim
!      fac2=fac1-fac
!      do 11 j=1,ndim
!        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
!11    continue

      if (ngtzero.eq.1) then
      ! only use above ptry if it is > 0. if < 0 calculate fac so that
      ! ptry = psum/(10*ndim)
      ! initialize
      facmin = fac
      facmax = fac
      small = 10.d0**(-30.d0)

       fac1=(1.-fac)/ndim
       fac2=fac1-fac
       do 11 j=1,ndim
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
         if (ptry(j).le.0.d0) then
          fact = (psum(j)*(1.d0-0.1d0)-p(ihi,j))
     &              /(psum(j)-(ndim+1)*p(ihi,j))
          if (fact.lt.facmin) facmin = fact
          if (fact.gt.facmax) facmax = fact
        end if
11    continue

! I don't know yet if facmin or facmax is good. try both if necessary
      fact = fac
      if (facmin.lt.fact) then
        write(*,*) 'trying facmin', facmin, fact
c       write(3,*) 'trying facmin', facmin, fact
        iwarn = 0
        fact = facmin
         fac1=(1.d0-fact)/ndim
         fac2=fac1-fact
        do j=1,ndim
           ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
          if (ptry(j).le.0.d0) iwarn = 1
        end do
        if (iwarn.eq.1) then
          write(*,*) 'facmin failed, trying facmax', facmax,fac
c         write(3,*) 'facmin failed, trying facmax', facmax,fac
          iwarn = 0
          fact = facmax
           fac1=(1.-fact)/ndim
           fac2=fac1-fact
          do j=1,ndim
             ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
            if (ptry(j).le.0.d0) then
               iwarn = 1
              ptry(j) = small
            end if
          end do
          if (iwarn.eq.1) write(*,*) 'facmax failed too... ptry->small'
c         if (iwarn.eq.1) write(3,*) 'facmax failed too... ptry->small'
        end if
      end if
      ! end more my addition
       else
       fac1=(1.d0-fac)/ndim
        fac2=fac1-fac
        do 17 j=1,ndim
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
17     continue
       end if

       ytry=funk(ptry)
       if (ytry.le.yb) then
         do 12 j=1,ndim
           pb(j)=ptry(j)
12      continue
         yb=ytry
       endif
       yflu=ytry-tt*dlog(ran1(idum))
       if (yflu.lt.yhi) then
         ! write(3,*) 'uphill accepted'
         y(ihi)=ytry
         yhi=yflu
         do 13 j=1,ndim
           psum(j)=psum(j)-p(ihi,j)+ptry(j)
           p(ihi,j)=ptry(j)
13      continue
       else
       ! write(3,*) 'uphill rejected'
       endif
       amotsa=yflu

       return
       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       FUNCTION ran1(idum)

c***********************************************************************
c      RAN1: Random number generator                                   *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - variable declarations changed from REAL to REAL*8             *
c                                                                      *
c      CM, Jan 2002                                                    *
c***********************************************************************

       INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
       REAL*8 ran1,AM,EPS,RNMX
       PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
!x     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
       INTEGER j,k,iv(NTAB),iy
       SAVE iv,iy
       DATA iv /NTAB*0/, iy /0/

       if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=NTAB+8,1,-1
           k=idum/IQ
           idum=IA*(idum-k*IQ)-IR*k
           if (idum.lt.0) idum=idum+IM
           if (j.le.NTAB) iv(j)=idum
11      continue
         iy=iv(1)
       endif

       k=idum/IQ
       idum=IA*(idum-k*IQ)-IR*k
       if (idum.lt.0) idum=idum+IM
       j=1+iy/NDIV
       iy=iv(j)
       iv(j)=idum
       ran1=min(AM*iy,RNMX)

       return
       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
