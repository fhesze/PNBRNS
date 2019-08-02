c    $Id: transcoeff-MT.f 16 2007-10-18 12:36:47Z centler $
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

       SUBROUTINE transcoeff_MT(k,aa,bb,c,di)

      include 'common_geo.inc'
       include 'common.inc'
       include 'common_drive.inc'

       real*8 DI(nx),AA(nx),C(nx),BB(nx),alphirr(nx)
      integer itransp,isolid
      real*8 vap,vam,va,ma,mb,mc,gradD,grade,gradv,ep,em,ea,gb,gac,
     &        gd, vvv, deltaxx, deltaxup, deltaxdown

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


c     alph0=200.D0                         !MT irrigation
      alph0=0.d0
      xirr=1/0.28                              !MT irrigation shallow
c     xirr=3.5
      do i=1,nx                              !MT irrigation
           alphirr(i)=alph0*exp(-x(i)/xirr)     !MT irrigation
c          alphirr(i)=0.
c          write(*,*) alphirr(i), i, x(i)
      enddo                                   !MT irrigation
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
        deltaxx=x(j+1)-x(j-1)
        deltaxdown=x(j)-x(j-2)
        deltaxup=x(j+2)-x(j)
        if (isolid.eq.1) then ! get the right advection velocity and eta
          eta=(1-por(j))
           vvv= vs(j)
           ai=0.
        else
           vvv= vd(j)
          eta=por(j)
           ai=alphirr(j)
        end if


        aa(j) = (-vvv/deltaxx - (disp(k,j)+disp(k,j-2))/2/
     &       (deltaxx*deltaxdown))*delt
        bb(j) = 1 + (vvv + ai*deltaxx + (disp(k,j)+disp(k,j-2))
     &               /2/deltaxdown
     &               +(disp(k,j)+disp(k,j+2))/2/deltaxup)*delt/deltaxx
         c(j) = -(disp(k,j)+disp(k,j+2))/2/(deltaxx*deltaxup)*delt
        di(j) = sp(k,j) + ai*spb(k,1)*delt

      end do

c***********************************************************************
c     boundary nodes                                                   *
c***********************************************************************
       do ii=1,2

        if (ii.eq.1) then
           j=1
           deltaxup=(x(3)-x(1))
        endif
        if (ii.eq.2) then
           j=nx
           deltaxdown=(x(nx)-x(nx-2))
        endif

c***********************************************************************
         ! fixed concentration
c***********************************************************************
         if (ibc(k,ii).eq.0) then
           aa(j) = 0.0d0
          c(j) = 0.0d0
          bb(j) = 1.0d0
          di(j) = spb(k,ii)
        else
              deltaxx=x(j+1)-x(j-1)

c***********************************************************************
           ! flux or gradient condition
c***********************************************************************
           if (isolid.eq.1) then     ! get the right velocity and eta
            vvv = vs(j)
           else
            vvv=vd(j)
           end if


c***********************************************************************
           ! known gradient
          if (ibc(k,ii).eq.1) then
             aa(j)= (-vvv/deltaxx - disp(k,j)/(deltaxx*2*(x(j)-x(j-1))))
     &               *delt
          bb(j)= 1 + (vvv/deltaxx + disp(k,j)/(deltaxx*2*(x(j)-x(j-1))))
     &               *delt
            c(j)=(-vvv/deltaxx-disp(k,j)/(deltaxx*2*(x(j+1)-x(j))))*delt
            di(j)=sp(k,j)+spb(k,ii)*disp(k,j)*2*delt/deltaxx
          end if
c***********************************************************************
          ! known flux (advective and diffusive)
          if(ibc(k,ii).eq.2) then
             if (j.eq.1) then
                aa(j)= 0.
                bb(j)= 1 + (vvv/deltaxx + (disp(k,j)+disp(k,j+2))/2
     &               /(deltaxx*deltaxup))*delt
               c(j) = -(disp(k,j)+disp(k,j+2))/2/(deltaxx*deltaxup)*delt
              di(j)=sp(k,j)+spb(k,ii)*delt/(eta*(x(j+1)-x(j-1)))
             else
               aa(j)=(-vvv/deltaxx - (disp(k,j)+disp(k,j-2))/2
     &                    /(deltaxx*deltaxdown))*delt
                bb(j)=1 + (disp(k,j)+disp(k,j-2))/2
     &                    /(deltaxx*deltaxdown)*delt
                c(j)=0.
                di(j)=sp(k,j)+spb(k,ii)*delt/(eta*(x(j+1)-x(j-1)))
             end if
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
