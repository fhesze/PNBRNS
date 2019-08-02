c    $Id: timestep.f 16 2007-10-18 12:36:47Z centler $
!        SUBROUTINE HAS BEEN TRANSLATED INTO FORTRAN 77
!    AND ADAPTED TO OUR VARIABLE DECLARATION  BY P. JOURABCHI
!  below: original header of subroutine received from Carl Steefel
!************** (C) COPYRIGHT 1995 Carl I. Steefel *******************
!
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:01
!
!                      All Rights Reserved
!
!  OSRT (Operator Splitting Reactive Transport) IS PROVIDED "AS IS"
!AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED. THE USER ASSUMES ALL RISKS
!  OF USING OSRT. THERE IS NO CLAIM OF THE MERCHANTABILITY OR FITNESS
!  FOR A PARTICULAR PURPOSE.
!
!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!  modified on 17.04.02
!**********************************************************************

      SUBROUTINE timestep()

      include 'common_geo.inc'
       include 'common.inc'
      include 'common_drive.inc'
      include 'timestep.inc'

      REAL*8 deltmax,theta,stepmax,dudtp,dudt,rnum,rden,dtemp
      INTEGER jx

      theta = 0.8d0
      deltmax = 2.0d0*delt !***MT
c     deltmax = 1.0d0*delt !***MT
!     stepmax = MIN(deltmax,dtmax,tstep)
      stepmax = MIN(deltmax,tstep)

      DO jx = 3,nx-2,2 !1,nx,2
           dudtp = (sp(kmast,jx) - spdt1(jx))/delt
           dudt = (spdt1(jx) - spdt2(jx))/dtold
           rnum = 2.0D0*theta*ttol*(spdt1(jx)+sp(kmast,jx))
           rden = (2.0D0/(delt+dtold))*(dudtp-dudt)
           IF ((dABS(rden) < 1.0d-30).or.
     &           (max(dabs(dudtp),dabs(dudt)).lt.1.d-10)) GO TO 20
          IF (((sp(kmast,jx)-spdt1(jx))/maxconc).lt.1.d-4) GO TO 20
           dtemp = DSQRT(DABS(rnum/rden))
           IF(dtemp/delt > 2.0D0) THEN
           dtemp = 2.0D0*delt
           END IF
           IF (dtemp*(1.0d0+1.0d-10) > stepmax) GO TO 15
           stepmax = dtemp
   15     CONTINUE
   20     CONTINUE
       END DO
      dtold = delt
      delt = stepmax

      RETURN
      END
!  *******************************************************

