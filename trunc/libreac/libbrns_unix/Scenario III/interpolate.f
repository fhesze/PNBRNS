c    $Id: interpolate.f 16 2007-10-18 12:36:47Z centler $

        subroutine interpolate(spint,j,t,k)

c***********************************************************************
c      INTERPOLATE:                                                    *
c      interpolate the calculated concentrations to the measurements   *
c      using splines.                                                  *
c                                                                      *
c      spint = interpolated concentration profile                      *
c      xx     = array with depth of calculated concentrations          *
c      xmeas = depth of measurements                                   *
c      in    = vector of calculated species (mapped from matrix sp)    *
c      y2    = second derivatives, used for spline                     *
c      spline: calculates splines (Press et al.)                       *
c      splint: uses splines to interpolate                             *
c      out   = interpolated concentration at depth xout                *
c                                                                      *
!!c    old comment:                                                    *
!!c      the interpolated concentrations are exported                  *
!!c      because splines don't work at the limits, it is checked if the*
!!c      measured and calculated concentrations are at the same depth  *
!!c      (more precisely within sx*deltax). if so, don't interpolate   *
!!c      sx    = fraction of dx within which no interpolation is made  *
!!c    check, I thought that should not really be a problem            *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************
c     uses spline and splint

         include 'common_geo.inc'
         include 'common.inc'
         include 'common_opt.inc'
         include 'common_meas.inc'
         dimension spint(maxxmeas),in(nx),y2(nx),xx(nx)
         real*8 spint,in,y2,xout,out
!!       real*8 sx
         integer i,j,t,k

!!        sx = 0.001d0
! space coordinate, but only at nodes (xx; x is at nodes and faces)
         do i=1,nx,2
           ii = (i+1)/2
!          x(ii) = (i-1) * delxi/2.d0
           xx(ii) = x(i)
         end do

         ! map "sp" onto "in"
         do i=1,nx,2
           ii = (i+1)/2
           in(ii) = sp(j,i)
         end do

! determine natural splines and evaluate the interpolated value -> Press
         nxtemp = int((nx-1)/2)
         call spline (xx,in,nxtemp,2.d30,2.d30,y2)
         do i=1,nrxmeas(t,k)
           xout = xmeas(t,k,i)
!!     ! check if measurement is at calculated depth (difference 
!! smaller than dx/1000
!!          ! if so pass spline. {floor(y) gives the nearest smaller 
!integer of y}
!!          if ( ((xout/delxi)-floor(xout/delxi)).ge.(sx*delxi)) then
!!            call splint (xx,in,y2,nxtemp,xout,out)
!!         ELSE
!!           OUT = IN(I)
!!          end if
           call splint (xx,in,y2,nxtemp,xout,out)
           spint(i) = out
         end do

         return
         end



! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       SUBROUTINE spline(xx,y,n,yp1,ypn,y2)

c***********************************************************************
c      SPLINE: Part of interpolation                                   *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - NMAX replaced by n, which is equal to ncomp                   *
c      - variable declarations changed from REAL to REAL*8             *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************

       INTEGER n,NMAX
       REAL*8 yp1,ypn,xx(n),y(n),y2(n)
       INTEGER i,k
       REAL*8 p,qn,sig,un,u(n)

       if (yp1.gt..99d30) then
         y2(1)=0.d0
         u(1)=0.d0
       else
         y2(1)=-0.5d0
         u(1)=(3.d0/(xx(2)-xx(1)))*((y(2)-y(1))/(xx(2)-xx(1))-yp1)
       endif

       do 11 i=2,n-1
         sig=(xx(i)-xx(i-1))/(xx(i+1)-xx(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(xx(i+
     *1)-xx(i))-(y(i)-y(i-1))/(xx(i)-xx(i-1)))/(xx(i+1)-xx(i-1))-sig*
     *u(i-1))/p
11    continue

       if (ypn.gt..99d30) then
         qn=0.
         un=0.
       else
         qn=0.5d0
         un=(3.d0/(xx(n)-xx(n-1)))*(ypn-(y(n)-y(n-1))/(xx(n)-xx(n-1)))
       endif
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

       do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
12    continue

       return
       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

       SUBROUTINE splint(xa,ya,y2a,n,xx,y)

c***********************************************************************
c      SPLINT: Part of interpolation                                   *
c      taken from Press et al. and documented there                    *
c      The following changes were made:                                *
c      - variable declarations changed from REAL to REAL*8             *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************

       INTEGER n
       REAL*8 xx,y,xa(n),y2a(n),ya(n)
       INTEGER k,khi,klo
       REAL*8 a,b,h

       klo=1
       khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.xx)then
           khi=k
         else
           klo=k
         endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.) pause 'bad xa input in splint'
       a=(xa(khi)-xx)/h
       b=(xx-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0

       return
       END

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
