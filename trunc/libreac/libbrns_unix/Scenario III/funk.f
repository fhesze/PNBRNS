c    $Id: funk.f 16 2007-10-18 12:36:47Z centler $
      function funk (ptry)

c***********************************************************************
c      FUNK:                                                           *
c      calculates the objective function associated with the vertex    *
c      supplied from amoeba (or amotry), called "ptry", which contains *
c      the values to the parameters being optimized                    *
c                                                                      *
c      FUNK contains calls to the main forward model and compares the  *
c      resulting calculated concentration profiles to the measured ones*
c      - starting time is 0 or the time of the last measurement to     *
c        which model calculations have just been compared              *
c      - the end time for simulation is the next time where            *
c        measurements are available                                    *
c      - parameters to be optimized are mapped on the entire paramter  *
c        list required in the forward code in the routines             *
c        "transferfw" and "transferback"                               *
c      - the calculated concentrations are interpolated to the depth   *
c        of the measured ones in the routine "interpolation"           *
c                                                                      *
c      CM, Dec 2001, modified Jan 2002                                 *
c***********************************************************************

      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

      real*8 ptry(nopt), spint(maxxmeas)
      real*8 of, funk, tstart, tend, av
      integer i,j,t
      real*8 tiny


c***********************************************************************
c      PARAMETER TRANSFER                                              *
c      put the optimized parameters into the complete parameter list   *
c     tiny    = small number, minimum value of any optimized parameter *
c               assumes that all parameters are positive               *
c               0 may lead to a division by 0 during optimization      *
c     ntotparam = total number of parameters                           *
c     nopt      = number of parameters which are optimized             *
c     idpar     = identification number of optimized parameters        *
c                 (i.e. index in entire parameter list)                *
c     ptry      = array of optimized parameters                        *
c     par       = array of all parameters                              *
c     transferback =  subroutine that associates par with the parameter*
c                     names in the maple/forward code                  *
c***********************************************************************
        tiny = 10.d0**(-30.)
           do k=1,ntotparam
             do i=1,nopt
             if (idpar(i).eq.k) then
              if (ngtzero.eq.1) then
              if(ptry(i).le.0.d0) then
                  write(3,*) 'set to tiny in funk'
               end if
              if(ptry(i).le.0.d0) ptry(i) = tiny ! force p > 0
             endif
             par(k) = ptry(i)
               end if
             end do
           end do
        call transferback() ! assign to their real maple names

c***********************************************************************
c      TIME LOOP                                                       *
c     tstart  = starting time of simulation. either 0 or last time     *
c               measurements have already been compared to calculated  *
c               values                                                 *
c     tend    = end time of next simulation. time of next measurement  *
c     ntopt   = number of different times measurements are available   *
c***********************************************************************
      of = 0.0d0 ! initialize objective function value

      do t=1,ntopt
        if (t.eq.1) then
          tstart = 0.0d0
        else
          tstart = timemeas(t-1)
        end if
        tend = timemeas(t)

c***********************************************************************
c      FORWARD CALCULATION                                             *
c***********************************************************************
        call diagenesis(tstart,tend)

c***********************************************************************
c      OBJECTIVE FUNCTION                                              *
c     of = value of the objective function                             *
c     nrspmeas(t) = number of species measured at a given time         *
c     spint = array of calculated concentrations interpolated to the   *
c            depths of measurements                                    *
c     funk = of                                                        *
c                                                                      *
c     write the parameters and the objective function to file no.3     *
c***********************************************************************
        do j=1,ncomp
          do k=1,nrspmeas(t)

            call objf(of,spint,j,k,t)

          end do ! end measured species loop
        end do   ! end component loop

      end do     ! end time loop

      funk = of


      write(*,13) of, (ptry(i),i=1,nopt)
      write(3,13) of, (ptry(i),i=1,nopt)
   13     format(100(e14.7,2x))

      end







