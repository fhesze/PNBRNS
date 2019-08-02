c    $Id: getdelt.f 16 2007-10-18 12:36:47Z centler $
      BLOCK DATA InitTimeStep
        common/tsparam/dtold,ttol,tstep,maxconc
        real*8 dtold,ttol,tstep,maxconc
         DATA dtold/1.d-6/, ttol/5.d-2/, tstep/1.0d-2/, maxconc/0.d0/
      END

      subroutine getdelt(nt,time,tend,spg)

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_drive.inc'
      include 'timestep.inc'

      integer jx,nt
      dimension spg(ncomp,nx)
      real*8 time,tend,dtt,dttest1,spg

c***********************************************************************
c  TIME STEP SELECTOR                                                  *
c     - the first time around the time step is either based on the     *
c     maple input or is selected based on diffusion/advection velocity *
c     i.e. the Courant-Friedrichs (CFL) and Nusselt (Nu) number        *
c     - subsequently, the timestep is either left at this initial value*
c     or selected based on an analysis of the second derivative of some*
c     "master" species, e.g. oxygen...                                 *
c                                                                      *
c     the initial timestep selection is based on the value of tfact:   *
c     tfact: if set to 0, dt defined in maple is used                  *
c            else, it is a safety margin and its value should be       *
c            between 0 and 1: the timestep based on CFL and Nu is      *
c            multiplied by tfact (delt=delt(CFL,Nu)*tfact              *
c     the subsequent timestepping is constant for nts=0, and depends   *
c     on the second derivative of the master species (indicated by the *
c     value of kmast) if nts is set to 1                               *
c                                                                      *
c     tfact, nts and kmast are all defined in drivervalue.f and are    *
c     transfered in common_drive.inc                                   *
c***********************************************************************

c***********************************************************************
c     INITIAL TIMESTEP                                                 *
c    based on Nu(=D*dt/(dx)^2) and CFL(=v*dt/dx) < 1 for stability,    *
c    where vd and D are total advection & diffusion coeff, respectively*
c    dttest,dttest1: temporary time step variables, initially a big nr.*
c    note that nt=0 before starting the timeloop, as soon as entering  *
c    it, nt=nt+1                                                       *
c                                                                      *
c    CM, Dec 2002                                                      *
c***********************************************************************
      if ((nt.eq.0).and.(tfact.gt.1.0d-12)) then
!     if ((nt.eq.0).and.(tfact.gt.0.0d0)) then

         dttest1 = 99999999.d0

        do j=2,nx-1,2
          if (dabs(vd(j)).gt.0.0d0) then
           dtt=(x(j+1)-x(j-1))/dabs(vd(j))
           if ((dtt.lt.dttest1).and.(dtt.gt.0.0d0)) dttest1=dtt
          end if
          do i=1,ncomp
             if (disp(i,j).gt.0.0d0) then
               dtt = (x(j+1)-x(j-1))*(x(j+1)-x(j-1))/disp(i,j)
               if ((dtt.lt.dttest1).and.(dtt.gt.0.0d0)) dttest1=dtt
            end if
           end do
        end do

         if (dttest1.gt.0.0d0) then
           dttest1 = dttest1 * tfact
           write(*,*) 'overwriting initial dt from...to:'
           write(*,*)  delt, dttest1
           delt = dttest1
         end if
      end if

c***********************************************************************
c     TIMESTEP ADJUSTMENT DURING SIMULATION                            *
c     based on looking at the change in the second derivative          *
c     first time around simply save the result for later calculation   *
c     of the derivatives. later during the simulation go into the      *
c     subroutine timestep to adjust timestep                           *
c                                                                      *
c     CM, July 2002 (work of Carl, original implementation here Parisa)*
c***********************************************************************
      if (nts.eq.1) then
        if (nt.le.1) then
          dtold = delt
          do jx = 1,nx,2
             spdt1(jx) = spg(kmast,jx)
          end do
        else
          do jx = 1,nx,2
            spdt2(jx) = spdt1(jx)
            spdt1(jx) = spg(kmast,jx)
            maxconc = max(sp(kmast,jx),maxconc)
          end do
           call timestep()
           delt=dmin1(dtold*1.001,delt) !*** MT
        end if
      end if

c***********************************************************************
c    AVOID JUMPING OVER FINAL TIME                                     *
c***********************************************************************
       if (tend.lt.(time+delt)) then
          write(*,*) 'restricting dt from...to...-> final time:'
          write(*,*)  delt, tend-time, tend
         delt=tend-time
       end if

      end

