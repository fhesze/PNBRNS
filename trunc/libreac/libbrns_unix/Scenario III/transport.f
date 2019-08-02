c    $Id: transport.f 16 2007-10-18 12:36:47Z centler $
c                  TRANSPORT
c***********************************************************************
c     species are transported unless the flag itransp is set .ne.0     *
c     the subroutine notransport is created from maple to indicate the *
c     species that are not transported                                 *
c     the transport is a linear problem, set up as                     *
c     Amat*Cvect = dvect, where Amat turns out to be tridiagonal       *
c     Amat and dvect are determined by the diffusion dispersion coeff.,*
c     the advection velocity, and the time step.                       *
c     see the subroutine transcoeff on how they are defined in detail  *
c***********************************************************************

#include<defines.inc>


       SUBROUTINE transport(k)

      include 'common_geo.inc'
       include 'common.inc'
       include 'common_drive.inc'

       real*8 DI(nx),AA(nx),C(nx),BB(nx)
      integer itransp

       itransp = 0            ! 0 indicates that species is transported
      call notransport(k,itransp)     ! set to 1 if not transported

      if (itransp.eq.0) then          ! regular transport
         call transcoeff_MT(k,aa,bb,c,di)     ! define Amat and d

#ifndef EXPLICIT_EULER
         call TRIDAG(AA,BB,C,di,co,nx)     ! solve Amat*C=d for C
#else
       write(*,*) "EXPLICIT!!"
       do j=3,nx-2,2                   ! explicit euler method
         co(j)=sp(k,j)-aa(j)*sp(k,j-2)-(bb(j)-1)*sp(k,j)-c(j)*sp(k,j+2)
       end do
         co(1)=sp(k,1) - (bb(1) - 1) * sp(k,1) - c(1) * sp(k,1+2)
         co(nx)=sp(k,nx) - aa(nx) * sp(k,nx-2) - (bb(nx) - 1) * sp(k,nx)
#endif
c override bottom boundary
c      co(nx)=co(nx-2)

      else                            ! no transport, but
        do j=1,nx,2                   ! transfer the current conc.
          co(j)=sp(k,j)
        end do
      end if

       RETURN
       END
