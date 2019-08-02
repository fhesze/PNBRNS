c    $Id: steadystate.f 16 2007-10-18 12:36:47Z centler $

       SUBROUTINE steadystate()

c***********************************************************************
c  calculation of steady state profiles                                *
c     - triggered by nsstate=1 or 2 (common.inc)                       *
c     0 = transport + reaction                                         *
c     system is solved for each component over space, and iterated     *
c      until all species are at steady state                           *
c     transport leads to matA(nx*nx)*Cnew(1*nx) = d_vec(1*nx)          *
c      - implemented either as Crank-Nicholson (see transcoeff.f)      *
c     reaction is implemented explicit, i.e. R = R(Cold)               *
c      this saves computer time, (reaction is only in d-vector)        *
c      but may lead to oscillating solutions. Therefore, I also use    *
c      the first derivative term dR/dC in the diagonal of matA         *
c      where R is the net rate acting on species C                     *
c      => transport  = reaction + d(reaction)/dC * (Cnew - Cold)       *
c      at steady state, Cnew = Cold, thus transport  = reaction        *
c      because many reaction kinetics acting on a species are driven   *
c      to a significant extent by the concentration of this species    *
c      (at least reactions consuming it), this seems to lead to a much *
c      more stable algorithm ("a bit implicit"). If higher order rate  *
c      laws are included, i.e. dR/dC = f(C), f(C) is f(Cold)           *
c                                                                      *
c      CM, March 02                                                    *
c***********************************************************************
      include 'common_geo.inc'
       include 'common.inc'
      include 'common_drive.inc'

      integer iflag(ncomp), iflagsum
      real*8 errmax, del, rat
       real*8 DI(nx),AA(nx),C(nx),BB(nx),ditot(nx)
      real*8 btot(nx), scal

         itss = 0
10        itss = itss+1          ! iteration index
          iflagsum = 0               ! # species not converged
         do k=1, ncomp
            itint = 0
18         itint = itint+1     ! internal loop index
           errmax = 0.d0 ! initialize maximum error found for species k
            do j=1,nx,2
              co(j) = sp(k,j)     ! initialize co,needed in transcoeff
            end do

           call transcoeff(k,aa,bb,c,di) ! matrix coeff. for transport

           ditot(1) = di(1)
           ditot(nx) = di(nx)
           btot(1) = bb(1)
           btot(nx)= bb(nx)

!          do j=3,nx-2,2                         ! preserve boundaries
           do j=1,nx,2
c***********************************************************************
c   SUBROUTINE ssrates                                                 *
c   rat  net rate acting on species k at node j, based on Cold         *
c   drdc is the derivative of the rate with respect to species k at j  *
c        calculated at Cold  = d(rat)/d(sp(k)) at j                    *
c   equilibrium reactions are simulated as kinetic ones,               *
c   using detailed balancing, i.e. Keq = [A]*[H+]/[HA]                 *
c   eq1: HA = A + Hplus: d[HA]/dt = -kf*[HA] + kb*[A]*[H] = 0          *
c   eq2: Keq = [A]*[H+]/[HA] (i.e. Keq is apparent eq. constant)       *
c   thus: kf/kb=[A]*[H]/[HA] = Keq => kb =kf/Keq                       *
c   a good choice for the value of kf is one that leads to kf*HA=1     *
c                                                                      *
c   CM, Mar 2002                                                       *
c***********************************************************************
              if ( ((j.eq.1).and.(ibc(k,1).eq.0))
     +          .or. ((j.eq.nx).and.(ibc(k,2).eq.0)) ) then
                !  preserve fixed conc @ boundaries
             else
c START ADDED
         isolid = 0
        call issolid(k,isolid)  ! 1 indicates solid
        if (isolid.eq.1) then! get the right advection velocity and eta
          ea = (1.0d0-por(j))   * area(j)
        else
          ea = por(j)   * area(j)
        end if
                call ssrates(rat,drdc,k,j)
               if (idrdc.eq.0) drdc = 0.d0
                ditot(j) = di(j)    - ea*rat  ! reaction term (explicit)
                ditot(j) = ditot(j) + ea*drdc*sp(k,j) ! include the
!                          derivate term
               btot(j)  = bb(j)    + ea*drdc
             end if
c END ADDED

c START REMOVED
c               call ssrates(rat,drdc,k,j)
c              if (idrdc.eq.0) drdc = 0.d0
c         ditot(j) = di(j)    - rat          ! reaction term (explicit)
c        ditot(j) = ditot(j) + drdc*sp(k,j) ! include the derivate term
c              btot(j)  = bb(j)    + drdc
c            end if
c END REMOVED
            end do
            call TRIDAG(AA,btot,C,ditot,co,nx) ! solve for new conc

           ! compute difference between old and new conc.
           scal = 0.d0
            do j=1,nx,2
             del = dabs(co(j)-sp(k,j))
             if(del.gt.errmax) errmax = del
              if (co(j).lt.0.d0) then          ! don't allow neg. conc.
                co(j)=0.d0
              end if
             scal = scal + (co(j)+sp(k,j))
            end do
           scal = scal/(dfloat(nx+1)/2.d0)/2.d0 ! results in average
!concentration

! debugging output (last and second last profile) if not converging
            if (itss.eq.itssmax) then
             do j=1,nx,2
               write(87,*) k, j, co(j), sp(k,j)
             end do
             write(87,*) ''
           end if

           ! update - only to a certain degree if new profile sucks
!           if ((errmax/scal.gt.0.001d0).and.(errmax.gt.etol)) then
!           if ((errmax/scal.gt.0.00001d0).or.(errmax.gt.etol)) then
            if ((errmax/scal.gt.0.00001d0).and.(errmax.gt.etol)) then
             iflag(k) = 1
              do j=1,nx,2
               if (itint.lt.itintmax) then
                 sp(k,j) = (ssrelax*co(j) + (1.d0-ssrelax)*sp(k,j))
               else
                 sp(k,j) = (co(j) + sp(k,j))*0.5d0
               end if
             end do
              if (itint.lt.itintmax) go to 18
           else
             iflag(k) = 0
           end if
           iflagsum = iflagsum + iflag(k)
          end do ! end species loop

         do k=1,ncomp
           do j=1,nx,2
             call limits(k,j)
           end do
         end do

c     ! enforce equilibrium
c         call eqrates()

         ! check if further iterations are needed
       if (((itss.lt.itssmax).and.(iflagsum.gt.0)).or.(itss.eq.1)) then
          goto 10
        end if

         if (itss.eq.itssmax) then
          write(*,*)'exceeding it. for steady state, not converged:'
c         write(5,*)'exceeding it., not converged species):'
          do k=1, ncomp
            write(*,*) k, iflag(k)            ! iflag=1: not converged
c            if (iflag(k).ne.0) write(5,*) k    ! iflag=1: not converged
          end do
        end if

        return
        end
