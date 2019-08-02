c    $Id: objf.f 16 2007-10-18 12:36:47Z centler $
      subroutine objf (of,spint,j,k,t)

c***********************************************************************
c      OBJF: calculates the objective function                         *
c                                                                      *
c      CM, Jan 2002                                                    *
c***********************************************************************
      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'

      real*8 spint(maxxmeas),of,std
      integer i,j,t,k

c***********************************************************************
c      OBJECTIVE FUNCTION                                              *
c     of = value of the objective function                             *
c    - defined as sum_times_measured(sum_space(deviations))            *
c      with the deviation defined as                                   *
c       {(measured-calculated conc.)                                   *
c        /average of the measurements over space at a given time}^2    *
c    - weighting the squared difference by the average scales the      *
c      different species or the standard deviation                     *
c    - the calculated concentrations are interpolated to the locations *
c      of measurements, using splines                                  *
c    nrspmeas(t) = number of species measured at a given time          *
c    idspmeas(t,k) = index of measured species k at a given time in    *
c                    the list of all calculated species                *
c    spint = array of calculated concentrations interpolated to the    *
c            depths of measurements                                    *
c    av = average measured concentration of a species at a given time  *
c    spmeas = measured concentration (time t, species k, depth i)      *
c    xmeas  = location of the measured concentration                   *
c    nrxmeas = number of measurements of the kth measured species      *
c             at time t                                                *
c                                                                      *
c    sigmeas: measured standard deviation, supplied in file            *
c             (getdat.f, not yet implemented)                          *
c    relsig: relative error e.g. 5% to estimate std.dev. if sigmeas is *
c            not supplied by the user                                  *
c    std: standard deviation used to calculate objective function      *
c    iof: toggle to choose between of(average) = 0 and of(sig) = 1     *
c                                                                      *
c    write the parameters and the objective function to file no.3      *
c***********************************************************************

          ! iof=0: scaling by average of measured conc.
         if (iof.eq.0) then
            if (idspmeas(t,k).eq.j) then ! was component measured ?
              call interpolate (spint,j,t,k) ! interpolate
              av = 0.d0               ! calculate the av. measured conc
             do i=1,nrxmeas(t,k)
                av = av+spmeas(t,k,i)
              end do
              av = av/nrxmeas(t,k)
              if (av.le.0.d0) av = 1.d0
              do i=1, nrxmeas(t,k) ! loop over all measured points
                of = of + ((spmeas(t,k,i)-spint(i))/av)**2. ! of
              end do
            end if
          end if

          ! iof=1: scaling by standard deviation of measured conc.
         if (iof.eq.1) then
            if (idspmeas(t,k).eq.j) then ! was component measured ?
              call interpolate (spint,j,t,k) ! interpolate
              do i=1, nrxmeas(t,k) ! loop over all measured points
               if (sigmeas(t,k,i).le.0.) then ! determine std.dev
                  std = relsig*spmeas(t,k,i)
               else
                 std = sigmeas(t,k,i)
               end if
               if (std.lt.stdmin) std = stdmin ! set minimum
                of = of + ((spmeas(t,k,i)-spint(i))/std)**2. ! of
             end do
            end if
          end if

      return
      end

