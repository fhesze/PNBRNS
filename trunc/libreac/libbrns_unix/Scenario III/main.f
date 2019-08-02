c    $Id: main.f 16 2007-10-18 12:36:47Z centler $
      program main
c***********************************************************************
c  this code is the fortran side of the biogeochemical reaction        *
c  network simulator, which consist of a maple interface and a fortran *
c  backbone for problem setup and solving, respectively                *
c  please do not alter the source code without talking to the authors  *
c  (see header below)                                                  *
c                                                                      *
c  credits for program parts taken or adapted from the literature are  *
c  given in the individual subroutines                                 *
c                                                                      *
c  July 2002, CM                                                       *
c***********************************************************************
      include 'common_geo.inc'
      include 'common.inc'
      include 'common_opt.inc'
      include 'common_meas.inc'
      include 'common_drive.inc'
c***********************************************************************
c  some explanations for variables defined in common blocks            *
c  - common_geo.inc                                                    *
c  nsolid:     number of solid species                                 *
c  ndiss:          number of dissolved species                         *
c  ncomp:          total number of species (nsolid+ndiss)              *
c  nreac:         total number of reactions                            *
c  nx:        gridcell number (note that concentrations are calculated *
c             at every second point, mixing parameters at the nodes    *
c             in between. node 1 and nx are the boundary points        *
c  kinetics:     kinetic rate constants and equilibrium constants      *
c  physics:     constants describing the physical environment          *
c             plus discretization                                      *
c  timestuff:     end time of simulation. value is defined in basic.f  *
c  - common.inc                                                        *
c  tspt:          sp, co and spold contain concentration profiles      *
c             disp, dsol_0, f_T deal with diffusion/dispersion coeff.  *
c               r are the reaction rates                               *
c  bound:     defines type of boundary conditions (value given in spb) *
c  - common_opt.inc                                                    *
c  nopt:          number of optimized parameters                       *
c  ntopt:          number of times measurements were made              *
c  ntotparam:     total number of parameters (see transferfw/back)     *
c  par:          entire list of parameters                             *
c  timemeas:     times measurements were taken                         *
c  idpar:          identifcation number of the optimized parameters    *
c  - common_meas.inc                                                   *
c  maxxmeas:     maximum numbers of measurements in a single profile   *
c  maxspmeas:     maximum number of species measured at the same time  *
c  idspmeas:     identification number of the species measured         *
c  nrxmeas:     number of measurements in a profile                    *
c  nrspmeas:     number of species measured at a given time            *
c  spmeas:     measured concentration profile                          *
c  xmeas:          depths of measurements                              *
c  sigmeas:     standard deviations of the measurements                *
c  - common_drive.inc                                                  *
c  see drivervalues.f                                                  *
c***********************************************************************

!!c      elapsed_time = TIMEF()

       write(*,*)
      write(*,*) '__________________ OPT-RTM-SED ______________________'
      write(*,*) '  reactive transport model with optimization layer   '
      write(*,*) '                   Version 1.0b                      '
      write(*,*) '_____________________________________________________'
      write(*,*) '  contributers:                                      '
      write(*,*) '  P.Regnier: basic RTM & concept                     '
      write(*,*) '  C.Meile: optimization, generalization&expansion RTM'
      write(*,*) '  D.Aguilera: web interface                          '
      write(*,*) '  P.Jourabchi: testing                               '
      write(*,*) '  All Rights Reserved                                '
      write(*,*) '_____________________________________________________'
      write(*,*)

       call printSvnVersion()

c***********************************************************************
c  DEFINITION OF DRIVER VARIABLES                                      *
c      defines the switches internal to fortran                        *
c***********************************************************************
       call drivervalues()

c***********************************************************************
c  START OF INITIALIZATION                                             *
c     sp        = concentration array  (derived from Maple 'variables')*
c     spb       = array of boundary cond. derived from Maple 'bnddata' *
c     basic     = routine defining Physical Parameters (Maple)         *
c     molecular     = routine defining the molecular diffusion coeffs. *
c                  and temperature dependance (Maple)                  *
c     boundaries = routine prescribing upper b.c.+initial cond.(Maple) *
c     biogeo     = routine defining the Biogeochemical param. (Maple)  *
c     initialcond= routine computing the initial profiles              *
c                  see routine drivervalues.f for the different options*
c                                                                      *
c     initial calls to get parameters. This is required if any of them *
c     are being optimized. if not this might be a waste of time ;-)    *
c     depth dependent profiles (porosity, area, D, v) are calculated   *
c     in diagenesis.f, because they may depend of the parameters to be *
c     optimized (e.g. dispersivity aL)                                 *
c***********************************************************************
c  obtain depth dependent parameters                                   *
c  this could be done in the main program only if porosity, area,      *
c  diffusion and advection velocities are not part of the optimization *
c  because I don't want to exclude this it is repeated here            *
c***********************************************************************
      call basic()
       call molecular()
       call biogeo()

       call gridsetup()
      call porarea()
      call advdiffcoeff()
      call printdepth()

       call boundaries()
       call initialcond()
cc     do jj=1,nx,2
cc 3011     format(200(1x,e14.7))
cc      write(232,3011) (sp(i,jj), i=1,ncomp)
cc      write(*,3011) (sp(i,jj), i=1,ncomp)
cc     pause
cc      end do
c***********************************************************************
c  START OF OPTIMIZATION                                               *
c     nopt = number of parameters to be optimized (Maple)              *
c            if nopt=0 bypass optimization and go for "final" run      *
c     transferfw: put all the parameters into an array called par      *
c     storedat: get the measured data the calculations are compared to *
c       -> getdata: filenames with data (Maple)                        *
c     optima: optimization of parameters using the downhill-simplex    *
c       -> funk: calculation of the objective function                 *
c           -> transferback: update parameters with the values in par  *
c           -> diagenesis: Pierre's core code (see below)              *
c           -> objf: calculates the objective function                 *
c              -> interpolation: interpolate calculated conc. to the   *
c                 measurements                                         *
c                  -> spline,splint: done with splines (Press et al.)  *
c       -> amoeba: simplex downhill optimization (see Press et al.)    *
c           -> amotry: part of amoeba                                  *
c              -> funk ...                                             *
c           -> funk     ...                                            *
c              -> transferback                                         *
c     optimlm: optimization of parameters using levenberg-marquardt    *
c           -> objf: -> interpolation: -> spline&splint                *
c     optftol: stopping criterium. details of use see amoeba & optimlm *
c                                                                      *
c     CM, Spring 2002                                                  *
c***********************************************************************
!!c     if (nopt.gt.0) then
      if (idoopt.eq.1) then

c***********************************************************************
c      TRANSFERFW: transfer maple & forward code parameters into a     *
c      generic list, called par                                        *
c      this list should contain all the parameters from common_geo.inc *
c      i.e. all rate parameters and all transport parameters           *
c      but not delt and delxi                                          *
c***********************************************************************
        call transferfw()

c***********************************************************************
c      STOREDAT: obtain the measured data                              *
c      read in the measured data in storedat.f                         *
c       the required format is described in the comments in storedat.f *
c      names of files containing measured data are listed in           *
c       the subroutine getdat.f, which is generated from within maple  *
c       and called from storedat                                       *
c                                                                      *
c      CM, Dec 2001                                                    *
c***********************************************************************
         call storedat()

c***********************************************************************
c      open file for output of the development of the objective        *
c      function and the parameter values, as they change along the way *
c      and do the optimization                                         *
c***********************************************************************
         open(unit=3,file='OF&par-values.dat',status='unknown')
        write(3,*) 'objective function & values of optimized parameters'
        if (iopt.eq.0) call optima()
        if (iopt.eq.1) call optimlm()
        if (iopt.eq.2) call optimsa()
         if (iopt.eq.3) call optimde()
        close(3)
      end if

c***********************************************************************
c  FINAL RUN                                                           *
c       diagenesis: forward model (Maple/Fortran)                      *
c       tstart: start time of simulation                               *
c       tend:     end time of simulation                               *
c***********************************************************************
      tstart = 0.d0
      tend = endt
       ntopt2 = 1
c     if ((ntopt.gt.0).and.(nopt.gt.0)) ntopt2 = ntopt
      if (idoopt.eq.1) ntopt2 = ntopt

c***********************************************************************
c  PRINT FINAL RUN                                                     *
c   limited to 100 components in format statement                      *
c***********************************************************************
       open(unit=5,file='conc.dat',status='unknown')
       open(unit=4,file='conc.txt',status='unknown')
! 3001     format(1x,f8.4, 1x, 200(1x,e14.7))
 3001     format(200(1x,e14.7))

      do ii=1,ntopt2
c        if ((ntopt.gt.0).and.(nopt.gt.0)) tend = timemeas(ii)
         if (idoopt.eq.1) tend = timemeas(ii)
         if ((tend.eq.0.).and.(nsstate.ne.1)) write(*,*) 'tend=0!'
        if (nsstate.eq.1) write(5,*) 'steady state, C(x), R(x)'
        if (nsstate.eq.2) write(5,*) 'steady state &transient C(x),R(x)'
        if ((nsstate.ne.1).and.(nsstate.ne.2)) then
           write(5,*) 'conc & rates @ time: ', tend
        end if
c       if (nsstate.eq.1) then
c          write(5,*) 'steady state concentration profile & rates'
c       else
c          write(5,*) 'conc & rates @ time: ', tend
c       end if
         call diagenesis(tstart,tend)
        do j=1,nx,2
           call rates(j)
           if(nsstate.eq.1) call rates(j)
          write(5,3001) x(j), (sp(i,j), i=1,ncomp),(r(k,j),k=1,nreac)
           ! create input file that can be read in initialcond.f
           if ((nsstate.eq.1).or.(nsstate.eq.2).or.(ii.eq.ntopt2)) then
c          if ((nsstate.eq.1).or.(ii.eq.ntopt2)) then
            write(4,3001) x(j), (sp(i,j), i=1,ncomp)
           end if
        end do
      end do
      close (4)
      close (5)

!!c      elapsed_time = TIMEF()
!!c      write(*,*) 'elapsed time: ', elapsed_time

      stop
      end






