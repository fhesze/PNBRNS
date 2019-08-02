c    $Id: optima.f 16 2007-10-18 12:36:47Z centler $
      subroutine optima()

c***********************************************************************
c      OPTIMA: Drives optimization using simplex downhill (AMOEBA)     *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

      integer iter, i, j, k, nrmeastot
      real*8 funk,p(nopt+1,nopt),psum(nopt),y(nopt+1)
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
c      choose the maximum fractional difference between best and       *
c      worst OF stored in simplex that leads to a stop of optimization *
c      optftol: some small number...                                   *
c      and call amoeba                                                 *
c***********************************************************************

        call amoeba(p,y,nopt+1,nopt,nopt,funk,iter)
        write(*,*) '# iterations performed:', iter


c***********************************************************************
c      put the optimized parameters back into the complete list        *
c      assign them to their real names (transferback)                  *
c      and give final output                                           *
c***********************************************************************
        do k=1,ntotparam
          do i=1,nopt
            if (idpar(i).eq.k) then
             par(k) = p(1,i)
              end if
          end do
        end do

        call transferback()


      write (*,*) 'BEST PARAMETER VALUES'
      write (3,*) 'BEST PARAMETER VALUES'

      do i=1,nopt
        write (*,*) i, par(idpar(i))
        write (3,*) i, par(idpar(i))
      end do

      return
      end

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       SUBROUTINE amoeba(p,y,mp,np,ndim,funk,iter)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

c***********************************************************************
c      AMOEBA: Parameter optimization using a simplex downhill method  *
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
c      CM, Dec 2001                                                    *
c***********************************************************************

        INTEGER iter,mp,ndim,np,t
        REAL*8 p(mp,np),y(mp),funk, fac
        EXTERNAL funk
CU     USES amotry,funk
        INTEGER i,ihi,ilo,inhi,j,m,n
        REAL*8 rtol,sum,swap,ysave,ytry,psum(np),amotry

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
         prob = 2.d0 !(any value greater than 1),so it doesn't interfere
       end if

        iter=0
 1     do 12 n=1,ndim
          sum=0.
         do 11 m=1,ndim+1
           sum=sum+p(m,n)
 11      continue
          psum(n)=sum
 12    continue
 2     ilo=1
        if (y(1).gt.y(2)) then
          ihi=1
          inhi=2
        else
          ihi=2
          inhi=1
        endif
        do 13 i=1,ndim+1
          if(y(i).le.y(ilo)) ilo=i
          if(y(i).gt.y(ihi)) then
            inhi=ihi
            ihi=i
          else if(y(i).gt.y(inhi)) then
            if(i.ne.ihi) inhi=i
          endif
 13    continue

        rtol=2.d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo)))
       if (iof.eq.1) prob = gammp(dfloat(idegf)/2.d0,y(ilo)/2.d0)
        if((rtol.lt.optftol).or.(prob.lt.problim)) then
!       if (rtol.lt.optftol) then
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
 14      continue
          return
        endif

        if (iter.ge.ITMAX) then
          write(*,*) 'ITMAX exceeded in amoeba'
          return
        end if

        iter=iter+2
        fac = -1.d0
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
        if (ytry.le.y(ilo)) then
          fac = 2.d0
          ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
        else if (ytry.ge.y(inhi)) then
          ysave=y(ihi)
          fac = 0.5d0
          ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
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
            iter=iter+ndim
            goto 1
          endif
        else
          iter=iter-1
        endif
        goto 2

        END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

c***********************************************************************
c      AMOTRY: Part of AMOEBA (simplex downhill)                       *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - common blocks added, so variables don't get undefined         *
c      - NMAX replaced by np, which is equal to ncomp                  *
c      - variable declarations changed from REAL to REAL*8             *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************

       INTEGER ihi,mp,ndim,np,j
       REAL*8 amotry,fac,p(mp,np),psum(np),y(mp),funk
       EXTERNAL funk
CU    USES funk
       REAL*8 fac1,fac2,ytry,ptry(np)
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
       fac1=(1.d0-fac)/ndim
       fac2=fac1-fac
       do 11 j=1,ndim
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
         if (ptry(j).le.0.d0) then
          fact = (psum(j)*(1.-0.1)-p(ihi,j))/(psum(j)-(ndim+1)*p(ihi,j))
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
          if (ptry(j).le.0.) iwarn = 1
        end do

        if (iwarn.eq.1) then
          write(*,*) 'facmin failed, trying facmax', facmax,fac
c         write(3,*) 'facmin failed, trying facmax', facmax,fac
          iwarn = 0
          fact = facmax
           fac1=(1.d0-fact)/ndim
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
       if (ytry.lt.y(ihi)) then
         y(ihi)=ytry
         do 12 j=1,ndim
           psum(j)=psum(j)-p(ihi,j)+ptry(j)
           p(ihi,j)=ptry(j)
12      continue
       endif
       amotry=ytry

       return
       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
