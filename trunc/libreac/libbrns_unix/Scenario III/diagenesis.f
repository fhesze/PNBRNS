c    $Id: diagenesis.f 16 2007-10-18 12:36:47Z centler $
#include<defines.inc>
       subroutine diagenesis (tstart,tend)

        include 'common_geo.inc'
        include 'common.inc'
        include 'common_drive.inc'

        dimension spguess(ncomp,nx)

      real*8 v_out, v_int
      real*8 tstart,tend,time
!     real*8 dtold
      integer nt

c***********************************************************************
c     FORWARD Reaction-Transport Code, based on operator splitting     *
c     and a reaction network created from within maple                 *
c     original version by Pierre                                       *
c     recent changes by CM:                                            *
c     - linked to optimization (interfacing)                           *
c     - included time step estimator based on CFL and Nu               *
c     - included dynamic timestepper (Parisa, Carl)                    *
c     - dealing differently with bad results                           *
c       - created warning flag                                         *
c       - set limits to catch NAN and out of bound concentration values*
c     - set depth to 0 at first node (not at dx/2)                     *
c     - changed sign of depth (i.e. distance from upper boundary)      *
c     - added depth dependent parameters (porosity, adv, disp.)        *
c     - made variable grid                                             *
c     - modified the transport matrix coefficients                     *
c       - to deal with variable grid                                   *
c       - defining fluxes at the boundary itself, not half a box inside*
c       - implemented different boundary conditions                    *
c     - also looping over boundary nodes                               *
c     - not reimposing boundary conditions after rxn                   *
c     - strang splitting, different solvers (strang, SNIA, steady st.) *
c     suggestions:                                                     *
c     + log transform concentrations in the reaction part              *
c                                                                      *
c     = check timestepper and interpolation                            *
c                                                                      *
c     CM, July 2002                                                    *
c***********************************************************************

c***********************************************************************
c  START OF INITIALIZATION                                             *
c  initalcond = routine prescribing the i.c. (maple)                   *
c     only used here if really starting at 0 again (i.e. time<delt)    *
c  obtain depth dependent parameters (porarea, advdiffcoeff)           *
c  this could be done in the main program only if porosity, area,      *
c  diffusion and advection velocities are not part of the optimization *
c  because I don't want to exclude this it is repeated here            *
c  however, this is only done when starting at t=0, because I assume   *
c  these properties not to vary with time (see governing equation)     *
c***********************************************************************
        if (tstart.lt.delt) then
         call porarea()
         call advdiffcoeff()
          call initialcond()
         call boundaries() ! added Jan'03
      ! set minimum conc to 1e-20
           do j=1,nx,2
                do i=1,ncomp
                     if(sp(i,j).lt.1.0d-20) sp(i,j)=1.0d-20
                enddo
           enddo
        end if

      do j=1,nx,2
           write(99,*) x(j),(sp(i,j),i=1,ncomp)
      enddo
c***********************************************************************
c  STEADY STATE - TRANSIENT SELECTOR                                   *
c    if nsstate = 1: steady state calculation                          *
c    if nsstate = 2: steady state calculation, then transient          *
c    else: transient only                                              *
c    nsstate is defined in drivervalues.f, not in maple                *
c    nssnow is set to 1 for ss calc (transcoeff)
c***********************************************************************
c       if (nsstate.eq.1) then
        if ((nsstate.eq.1).or.(nsstate.eq.2)) then
         nssnow=1 ! steady state calculation (for transcoeff.f)

        if(istst.eq.1)call steadystate() ! transp(C-N) = Rxn + drdc*delC
!        if(istst.eq.2)call globess()
!        if(istst.eq.3)call steadystate2() ! jacobian, func, all in one

         if (nsstate.eq.1) return
c        return
       end if
       nssnow = 0 ! indicates transient calculation

c***********************************************************************
c  initial TIME STEP SELECTOR                                          *
c  deltsave saves the original (maple) timestep, which is restored     *
c  before leaving this subroutine. needed if delt is changed during    *
c  the simulation (e.g. to avoid jumping over tend)                    *
c***********************************************************************
      deltsave = delt
      nt = 0
      time = tstart
      call getdelt(nt,time,tend,spguess)
!     call getdelt(nt,time,tend,dtold,spguess)

c***********************************************************************
c  START OF TIME LOOP                                                  *
c***********************************************************************
      do while (time.lt.tend)
999       nt = nt+1
        time = time + delt
c       write(*,*) time, delt
       if (time.le.v_out.and.v_out.lt.time+delt) write(*,*) time

c***********************************************************************
c  START OF TRANSPORT                                                  *
c     ncomp     = total number of variables (Maple)                    *
c     nx        = total number of grid points (node&faces) (maple)     *
c     co        = concentration array (temporary array)                *
c     spguess   = concentration prior to the transport step            *
c                 (used in the batch reactor procedure)                *
c     transport = routine for advection/ dispersion scheme             *
c***********************************************************************
         do k = 1,ncomp
           do j = 1,nx,2
             co(j) = sp(k,j)
             spguess(k,j) = sp(k,j)
           end do
           call transport(k)
           do j = 1,nx,2
             sp(k,j) = co(j)
           end do
c          pause 'end of transport'
         end do
c       pause
c***********************************************************************
c  START OF REACTION                                                   *
c     fmax  = maximum calculated residual in current time step         *
c             (residual x time step) for all components in the entire  *
c             spatial domain                                           *
c     emax  = maximum calculated change in concentration for all       *
c             components in the entire spatial domain.                 *
c     spold = the concentration profile immediately after transport    *
c     iflag = flag to indicate problems in solver                      *
c             (0=no problem, 1=negative concentrations,                *
c              2=exceeding iterations, 3= 1&2)                         *
c             warning is given only once for all species and all nodes *
c     isol  = selector for solver                                      *
c             (1=newton with relaxation, 2=newton with linesearch)     *
c     newtonsub = newton root finding (isol = 1)                       *
c     newt      = newton with additional linesearch (isol = 2)         *
c***********************************************************************
         iflag = 0
         do k=1,ncomp
           do j = 1,nx,2
            spold(k,j) = sp(k,j)
            sp(k,j) = spguess(k,j)
          end do
         end do

         depth =  0.d0 ! depth in the sediment column
!       call rates(1)
!        call out(1,nt,time,depth,v_out,v_int)

!       do j = 3,nx-2,2  ! outer space loop for reaction network
        newtonflag=0          !MT
        newtoncounter=0     !MT
         do j = 1,nx,2    ! outer space loop for reaction network
c          write(*,*) j
         if (isol.eq.1) then
           call newtonsub(j,iflag,ftol,etol,newton,ilud)
         end if
         if (isol.eq.2) call newt(check,j,iflag,ftol,etol,newton,ilud)
           if (iflag.eq.2) then                    !MT
                newtonflag=2                         !MT
                newtoncounter=newtoncounter+1     !MT
                write(*,*) "Exceeding newton iterations in node:"
                write(*,*) j
           endif                                        !MT
c***********************************************************************
c  DEAL WITH BAD RESULTS                                               *
c     if a concentration gets VERY large or negative, set concentration*
c     back to ±1e30. This can happen during optimization if the        *
c     parameters are chosen very poorly. clumsy if statements should   *
c     also set the concentration to these values if one runs into a NAN*
c***********************************************************************
           do i=1,ncomp
            call limits(i,j)
!            if ((sp(i,j).lt.1.0d30).and.(sp(i,j).gt.-1.0d30)) then
!           ! do nothing. but if NAN enter in ...else...
!            else
!              write(*,*) 'out of bounds...reset to +/-1e30'
!              if(sp(i,j).gt.0.) then
!                sp(i,j) = 1.0d30
!              else
!                if (sp(i,j).lt.0.d0) then
!                  sp(i,j) = -1.0d30
!                else
!                  sp(i,j) = 0.d0
!                end if
!              end if
!            end if
           end do

         end do  ! end of space loop

         if (iflag.eq.1) write(*,*) 'negative concentrations!'
         if (newtonflag.eq.2) then          !MT
           write(*,*) 'exceeding newton iterations'
           write(*,*) time, delt, newtoncounter     ! MT
           time=time-delt               ! MT
           nt=nt-1                         ! MT
           delt=delt/2               ! MT
           goto 999                    ! MT
        endif
         if (iflag.eq.3) write(*,*) 'conc<0 & exceeding newton it.'

c***********************************************************************
c   transport again if doing strang splitting                          *
c   - update co with the current concentrations                        *
c   - send to transport (timestep is adapted in transcoeff.f)          *
c   - update concentrations                                            *
c***********************************************************************
         if (isplit.eq.1) then
           do k = 1,ncomp
             do j = 1,nx,2
               co(j) = sp(k,j)
             end do
             call transport(k)
             do j = 1,nx,2
               sp(k,j) = co(j)
             end do
           end do
        end if

c***********************************************************************
c   OUTPUT                                                             *
c     v_out  = first time for writing to data file                     *
c     v_int  = time interval between data output                       *
c     out    = routine handling the output files (Maple)               *
c***********************************************************************
          do j = 1,nx,2
            call rates(j)
            call out(j,nt,time,x(j),v_out,v_int)
         end do

c***********************************************************************
c   update timestep                                                    *
c***********************************************************************
        if (time.lt.tend) then
          call getdelt(nt,time,tend,spguess)
!         call getdelt(nt,time,tend,dtold,spguess)
        end if

       end do          ! end of time loop
      delt = deltsave ! restore the original (maple-) timestep

       return
       end
