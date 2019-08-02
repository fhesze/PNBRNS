c    $Id: drivervalues.f 16 2007-10-18 12:36:47Z centler $
        subroutine drivervalues()
c***********************************************************************
c     this subroutine pools (most of) the "hidden" variables that      *
c     define the different options in the code,                        *
c     but which are not explicitly defined in the maple interface      *
c                                                                      *
c     CM, Spring 2002                                                  *
c***********************************************************************
         include 'common_geo.inc'
         include 'common_opt.inc'
         include 'common_drive.inc'

c***********************************************************************
c      steady state or transient forward simulation                    *
c        nsstate = 1:  steady state calculation                        *
c        nsstate = 2:  steady state calculation followed by transient  *
c        other value for nsstate: transient simulation                 *
c        uc     spatial discretization scheme:        *
c                          = 0 for unknown scheme, otherwise    upwind *
c        alphaCN:      implicit explicit weighting. must be >=0, <=1   *
c                      0.5 is Crank-Nicholson, 1 is implicit transport *
c                      if nsstate =1, you might want alphaCN=1.0d0     *
c                      ...alphaCN=0 leads to BET=0 in steadystate mode *
c                         in such a case, alphaCN is set to 1          *
c        isplit        type of operator splitting                      *
c                      default (i0e. <> 1) sequential non-iterative,   *
c                       i.e. transport, then reaction                 *
c                      1: strang (half transp., rxn, half transport)   *
c                       strang might be better for enforcing b. cond.  *
c                       and reduce splitting error but costs more doing*
c                       the transport                                  *
c         ilud:     1: LU-decomposition following Press                *
c                    else: LU-decomposition using lapack               *
c***********************************************************************
          nsstate = 0 ! 2
          uc = 0.5d0 ! 0.0d0
         alphaCN = 1.0d0 !0.5d0 ! 1.0d0
         isplit = 0
         ilud = 1

c     if ((nsstate.eq.1).and.(alphaCN.eq.0.0d0)) then
      if (((nsstate.eq.1).or.(nsstate.eq.2)).and.(alphaCN.eq.0.d0)) then
        write(*,*) 'steady state and alphaCN=0 is incompatible'
        write(*,*) 'alphaCN set to 1'
         alphaCN = 1.0d0
      end if

c***********************************************************************
c      iteration for steady state (RELEVANT FOR STEADY STATE ONLY)     *
c        itssmax:          maximum number of iteration loops (~5000)   *
c        itintmax:     number of loops per species before going to next*
c                     iteration (internal loops due to explicit term   *
c                     in rate description) (~10)                       *
c        idrdc:          if 1 then use d(rate)/d(C) term to calculate  *
c                     steady state. if 0 then don't. 1 is recommended, *
c                     for details see steadystate.f                    *
c        ssrelax:          relaxation term (amount of suggested step to*
c                     improve the solution you actually take) (~0.3)   *
c        istst:         1: steadystate.f (iterative by species)        *
c                     2: globe.f (global assembly, full dR/dC)         *
c                         current version only allows for istst=1      *
c with Parisa's global addition, one will have to change istst here,   *
c add the subroutines globem and globess and uncomment the call to     *
c globe.f in diagenesis.f                                              *
c also, the processor file will need to be modified, uncommenting      *
c the call to acg18(rjac,dir_f)                                        *
c***********************************************************************
       itssmax = 5000 ! 5000
       itintmax = 1 ! 10
       idrdc = 1
       ssrelax = 1.0d0 ! 0.3d0
       istst=1

       if(istst.ne.1)then
           write(*,*) 'istst must be 1: drivervalues.f'
           stop
       endif
!      if((istst.lt.1).or.(istst.gt.2))then
!          write(*,*) 'istst must be 1 or 2: drivervalues.f'
!         pause
!       end if

c***********************************************************************
c      time stepping parameter (getdelt.f) (TRANSIENT)                 *
c        tflag = 0.0d0: use the timestep set in the maple input sheet  *
c        tflag > 0: delta_t is based on tflag*(CFL,Nu) criterion       *
c        nts = 1: update timestep during simulation based on second    *
c                 derivative of the master species (see kmast)         *
c        nts = 0: don't update timestep                                *
c        kmast: indicates master species whose second derivative with  *
c               respect to time is used to adjust timestep             *
c***********************************************************************
          tfact = 0.0d0
         nts = 0
         kmast  = 1

         if ((tfact.gt.1.d0).or.(tfact.lt.0.d0)) then
           write(*,*) 'tfact must be between 0 and 1: drivervalues.f'
           pause
         end if
         if ((nts.eq.0).or.(nts.eq.1)) then
            if ((kmast.lt.1).or.(kmast.gt.ncomp)) then
             write(*,*) 'kmast must be between 1&ncomp: drivervalues.f'
             pause
            end if
         else
           write(*,*) 'nts must be 0 or 1: drivervalues.f'
           pause
         end if

c***********************************************************************
c      solver parameters (TRANSIENT)                                   *
c      - choice of solver in diagenesis.f                              *
c        isol=1: newton with relaxation (newtonsub.f)                  *
c        isol=2: newton with linesearch (newt.f)                       *
c      - newton: number of newton iterations                           *
c      - etol:   absolute difference allowed in convergence            *
c      - ftol:   relative difference allowed in convergence            *
c        (for details on etol and ftol see newtonsub.f and newt.f)     *
c***********************************************************************
          isol = 1
          newton = 100000 !*10 ! 100
          etol = 1.d-12 ! 6
         ftol = 5.d-12

      if ((isol.lt.1).or.(isol.gt.2)) then
        write(*,*) 'isol set to 1'
        isol = 1
      end if

c***********************************************************************
c      optimization parameters                                         *
c      iopt: selects optimization routine                              *
c            0: simplex downhill                                       *
c            1: levenberg-marquardt                                    *
c            2: simulated annealing                                    *
c            3: differential evolution (genetic algorithm family)      *
c      optftol: defines optimization stop criterion: delt_of/magnitude *
c      problim: probability to accept a wrong model (stop criterium)   *
c      perturb: perturbation of vertices: p = p*(1+perturb)            *
c               for LM (iopt=1), this should be rather small (0.01)    *
c               for DS & SA (iopt=0,2), 0.1 seems ok                   *
c      itmax: maximum number of iterations in optimization             *
c      iof: selects type of objective function (objf.f)                *
c           0: weighting by average of measured profile                *
c           1: weighting by uncertainty (std.dev)                      *
c      relsig: relative error: std.dev = relsig*Cmeas                  *
c      stdmin: minimum error in Cmeas                                  *
c      ngtzero: if 1, then parameters are forced to be > 0             *
c                                                                      *
c      NOTE: there are some local switches in                          *
c       optimsa (cooling schedule: ncool, alphatemp)                   *
c       optimde (crossover probabilities, range, a.s.o)                *
c***********************************************************************
        iopt = 3
        optftol =10.d0**(-9.d0)
        problim = 0.001d0
       perturb = 0.1d0
       itmax = nopt*50
       iof = 1
       relsig = 0.01d0
       stdmin = 10.0d0**(-6.d0)
       ngtzero=0

       ! obj.function must be related to std.dev in LM
       if (iopt.eq.1) iof = 1

        return
       end
