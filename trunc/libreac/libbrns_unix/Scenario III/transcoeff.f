c    $Id: transcoeff.f 16 2007-10-18 12:36:47Z centler $
c                  TRANSPORT coefficients

c***********************************************************************
c     the transport equation is discretized and written in matrix form *
c     Amat*Cvect = dvect, where Amat turns out to be tridiagonal       *
c     the steady state coefficients are similar to the transient ones  *
c     except that components due to the time derivative are not added  *
c     see separate documentation file for details what is solved (&how)*
c                                                                      *
c     CM, July 2002                                                    *
c***********************************************************************

       SUBROUTINE transcoeff(k,aa,bb,c,di)

      include 'common_geo.inc'
       include 'common.inc'
       include 'common_drive.inc'

       real*8 DI(nx),AA(nx),C(nx),BB(nx)
      integer itransp,isolid
      real*8 vap,vam,va,ma,mb,mc,gradD,grade,gradv,ep,em,ea,gb,gac,gd

c***********************************************************************
c     no transport ?                                                   *
c     if itransp is <> 0 after the call to notransport, the species is *
c     not being transported at all. If so, set all coefficients to 0,  *
c     note that in the transient code you never get here for a species *
c     that is not transported (saves time), and it is done here only   *
c     for the steady state problem                                     *
c                                                                      *
c     CM, July 2002                                                    *
c***********************************************************************
       itransp = 0
      call notransport(k,itransp)  ! itransp set to 1 if not transported
      if (itransp.ne.0) then
!        do j=1,nx,2
         do j=3,nx-2,2       ! only interior, else tridag fails
          aa(j) = 0.0D0
           bb(j) = 1.0d0 ! 0.0D0
           c(j)  = 0.0D0
           di(j) = co(j) ! 0.0d0
!          conserve initial cond. at boundary
!          if ((j.eq.1).or.(j.eq.nx)) then
!            bb(j) = 1.0d0
!           di(j) = sp(k,j)
!          end if
         end do
        bb(1) = 1.0d0
        di(1) = sp(k,1)
        aa(1) = 0.d0
        c(1)  = 0.d0
        bb(nx) = 1.0d0
        di(nx) = sp(k,nx)
        aa(nx) = 0.d0
        c(nx)  = 0.d0
        return
       end if

c***********************************************************************
c     regular transport                                                *
c     - half timestep if doing strang splitting                        *
c     - check if a species is a solid (isolid=1) or a solute,          *
c       because this has an effect on the selection which transport    *
c       parameters and theta (por*area or (1-por)*area) is used        *
c     - then calculate the coefficients in the interior of the domain  *
c     - then calculate the boundary nodes                              *
c     - add time dependence                                            *
c     - double timestep again (restore) if doing strang splitting      *
c                                                                      *
c     CM, July 2002                                                    *
c***********************************************************************


c***********************************************************************
c     half time step if doing strang splitting (restore at the bottom) *
c***********************************************************************
       if (isplit.eq.1) then
        delt = delt*0.5d0
      end if

c***********************************************************************
c     solid species?                                                   *
c***********************************************************************
       isolid = 0
      call issolid(k,isolid)      ! 1 indicates solid

c***********************************************************************
c     interior nodes                                                   *
c***********************************************************************
       do j=3,nx-2
        if (isolid.eq.1) then ! get the right advection velocity and eta
          vap = vs(j+1)
          vam = vs(j-1)
          ep = (1.0d0-por(j+1)) * area(j+1)
          em = (1.0d0-por(j-1)) * area(j-1)
          ea = (1.0d0-por(j))   * area(j)
        else
           vap = vd(j+1)
           vam = vd(j-1)
          ep = por(j+1) * area(j+1)
          em = por(j-1) * area(j-1)
          ea = por(j)   * area(j)
        end if

      if (uc.ne.0.d0) then
c     ***************** simple upwind scheme ****************

        ma =  disp(k,j-1)*em/(x(j)-x(j-2)) /(x(j+1)-x(j-1))
     +         + vam*em/(x(j)-x(j-2))
     +
        mb = ( -disp(k,j+1)*ep/(x(j+2)-x(j))
     +         -disp(k,j-1)*em/(x(j)-x(j-2)) )/(x(j+1)-x(j-1))
     +         -((vap+vam)/2.0d0)*ea/(x(j)-x(j-2))

        mc = ( disp(k,j+1)*ep/(x(j+2)-x(j)) )/(x(j+1)-x(j-1))
      else

c      unknown scheme which involved estimating fluxes at faces
c      and is described in TM2.doc

        ma = ( disp(k,j-1)*em/(x(j)-x(j-2))
     +         + vam*em*(x(j)-x(j-1))/(x(j)-x(j-2)) )
     +       /(x(j+1)-x(j-1))
        mb = ( -disp(k,j+1)*ep/(x(j+2)-x(j))
     +         -disp(k,j-1)*em/(x(j)-x(j-2))
     +         -vap*ep*(1.0d0 - (x(j+1)-x(j))/(x(j+2)-x(j)))
     +         +vam*em*(1.0d0 - (x(j)-x(j-1))/(x(j)-x(j-2))) )
     +       /(x(j+1)-x(j-1))
        mc = ( disp(k,j+1)*ep/(x(j+2)-x(j))
     +         - vap*ep*(x(j+1)-x(j))/(x(j+2)-x(j)) )
     +       /(x(j+1)-x(j-1))
!     +         + vap*ep*(x(j+1)-x(j))/(x(j+2)-x(j)) )

      end if

         aa(j) = alphaCN*ma
        bb(j) = alphaCN*mb
        c(j)  = alphaCN*mc
        di(j) = (alphaCN - 1.0d0) * (ma*co(j-2) + mb*co(j) + mc*co(j+2))

c***********************************************************************
        ! transient
c***********************************************************************
c       if (nsstate.ne.1) then
        if (nssnow.eq.0) then
c         bb(j) = bb(j) - por(j)*area(j)/delt
c         di(j) = di(j) - por(j)*area(j)/delt*co(j)
          bb(j) = bb(j) - ea/delt
          di(j) = di(j) - ea/delt*co(j)
         end if
      end do

c***********************************************************************
c     boundary nodes                                                   *
c***********************************************************************
       do ii=1,2

        if (ii.eq.1) j=1
        if (ii.eq.2) j=nx

c***********************************************************************
         ! fixed concentration
c***********************************************************************
         if (ibc(k,ii).eq.0) then
           aa(j) = 0.0d0
          c(j) = 0.0d0
          bb(j) = 1.0d0
          di(j) = spb(k,ii)
        else

c***********************************************************************
           ! flux or gradient condition
c***********************************************************************
           if (isolid.eq.1) then     ! get the right velocity and eta
            vap = vs(j+1)
            vam = vs(j-1)
            va = vs(j)
            ep = (1.0d0-por(j+1)) * area(j+1)
            em = (1.0d0-por(j-1)) * area(j-1)
            ea = (1.0d0-por(j))   * area(j)
           else
             vap = vd(j+1)
             vam = vd(j-1)
             va = vd(j)
            ep = por(j+1) * area(j+1)
            em = por(j-1) * area(j-1)
            ea = por(j)   * area(j)
           end if
          ! calculate the gradients of D, eta and v at the boundary
          gradD = ( -(x(j+1)-x(j)) * disp(k,j-1)
     +           +(x(j+1)-x(j)-x(j)+x(j-1)) * disp(k,j)
     +           +(x(j)-x(j-1)) * disp(k,j+1) )
     +          /( (x(j)-x(j-1))*(x(j+1)-x(j-1)) )
          grade = ( -(x(j+1)-x(j)) * em
     +           +(x(j+1)-x(j)-x(j)+x(j-1)) * ea
     +           +(x(j)-x(j-1)) * ep )
     +          /( (x(j)-x(j-1))*(x(j+1)-x(j-1)) )
          gradv = ( -(x(j+1)-x(j)) * vam
     +           +(x(j+1)-x(j)-x(j)+x(j-1)) * va
     +           +(x(j)-x(j-1)) * vap )
     +          /( (x(j)-x(j-1))*(x(j+1)-x(j-1)) )

!    calculate matrix coefficients assuming the ghost points are real
!    because ghost point is at j+1 or j-1, not j+2 or j-2,
!    (i.e. 0 and nx+1 are ghost points, which are both node and face)
!    rewrite ma,mb,mc
           if (j.eq.1) then
            xp = x(j+2)
            xm = x(j-1)
          else
            xp = x(j+1)
            xm = x(j-2)
          end if

          ma = ( 2.0d0*disp(k,j)*ea
     +           + (disp(k,j)*grade +ea*gradD -va*ea)*(-xp+x(j)) )
     +         /((x(j)-xm)*(xp-xm))
          mb = ( -2.0d0*disp(k,j)*ea
     +           + (disp(k,j)*grade +ea*gradD -va*ea)
     +              *(xp-x(j)-x(j)+xm) )
     +         /((x(j)-xm)*(xp-x(j)))
     +         + (-va*grade -ea*gradv)
          mc = ( 2.0d0*disp(k,j)*ea
     +           + (disp(k,j)*grade +ea*gradD -va*ea)*(x(j)-xm) )
     +         /((xp-x(j))*(xp-xm))


!    calculate the contribution to the different nodes
!    when expressing the ghost point with known quantities
           if (ii.eq.1) then
            xout = x(0)
            xin  = x(3)
          else
            xout = x(nx+1)
            xin  = x(nx-2)
          end if

c***********************************************************************
           ! known gradient
          if (ibc(k,ii).eq.1) then
            gac = (xout-x(j))*(xout-x(j))/((x(j)-xin)*(x(j)-xin))
            gb =-(xout-x(j)-x(j)+xin)*(xout-xin)/((x(j)-xin)*(x(j)-xin))
            gd  = spb(k,ii)*(xout-x(j))*(xout-xin)/(x(j)-xin)
          end if
c***********************************************************************
          ! known flux (advective and diffusive)
          if(ibc(k,ii).eq.2) then
            gac = (xout-x(j))*(xout-x(j))/((x(j)-xin)*(x(j)-xin))
            gb  = ( -disp(k,j)*ea*(xout-x(j)-x(j)+xin)
     +               /((x(j)-xin)*(xout-x(j)))
     +              - disp(k,j)*grade - ea*gradD + va*ea )
     +            * (xout-x(j))*(xout-xin)/(disp(k,j)*ea*(x(j)-xin))
            gd  = -spb(k,ii)*(xout-x(j))*(xout-xin)
     +            /(disp(k,j)*ea*(x(j)-xin))
          end if

           if (j.eq.1) then
            cout = gb*co(j) + gac*co(j+2) + gd
            aa(j) = 0.0d0
            bb(j) = alphaCN * (mb + ma*gb)
            c(j)  = alphaCN * (mc + ma*gac)
            di(j) = (alphaCN-1.0d0) * (ma*cout + mb*co(j) + mc*co(j+2))
            di(j) = di(j) - alphaCN*ma*gd
          else
            cout = gb*co(j) + gac*co(j-2) + gd
            aa(j) = alphaCN * (ma + mc*gac)
            bb(j) = alphaCN * (mb + mc*gb)
            c(j)  = 0.0d0
            di(j) = (alphaCN-1.0d0) * (ma*co(j-2) + mb*co(j) + mc*cout)
            di(j) = di(j) - alphaCN*mc*gd
          end if

           if ((uc.ne.0.d0).and.(ibc(k,ii).eq.1)) then
            ! if (j.eq.1) then
            !  only for the lower boundary (often a zero gradient)
            ! end if
            if (j.eq.nx) then
              aa(j) = alphaCN*(1.d0)
              bb(j) = alphaCN*(-1.d0)
              di(j) = (alphaCN-1.0d0) * (co(j-2) - co(j))
              di(j) = di(j) - spb(k,ii)*(x(j)-x(j-2))
            end if
          end if


      !! modification for gradient condition and steady state solver
      !! only do it for gradient, not flux, because flux at bottom of 
      !! solids is advection dominated - problem seems diffusion related
           if ((nssnow.ne.0).and.(ibc(k,ii).eq.1)) then
            ! if (j.eq.1) then
            !!leave at this time, because often the upper boundary
            !!involves a lot of reaction, hence gradient estimate across
            !!the boundary is a good thing
            ! end if
            if (j.eq.nx) then
              aa(j) = alphaCN*(-1.d0)
              bb(j) = alphaCN*(1.d0)
              di(j) = (alphaCN-1.0d0) * (-co(j-2) + co(j))
              di(j) = di(j) + spb(k,ii)*(x(j)-x(j-2))
            end if
          end if


c***********************************************************************
          ! transient (and flux or gradient condition)
c***********************************************************************
c          if (nsstate.ne.1) then
           if (nssnow.eq.0.and.uc.eq.0.d0) then
c           bb(j) = bb(j) - por(j)*area(j)/delt
c           di(j) = di(j) - por(j)*area(j)/delt*co(j)
            bb(j) = bb(j) - ea/delt
            di(j) = di(j) - ea/delt*co(j)
           end if

        end if
      end do

c***********************************************************************
c     restore time step if doing strang splitting                      *
c***********************************************************************
       if (isplit.eq.1) then
        delt = delt*2.0d0
      end if

!     if (k.eq.1) then
!       write(*,*) 'show'
!     if (co(nx-2).lt.0.505d0) then
!            write(*,*) 'show'
!      end if
!     end if


       RETURN
       END
