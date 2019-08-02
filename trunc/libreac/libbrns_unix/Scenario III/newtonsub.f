c    $Id: newtonsub.f 16 2007-10-18 12:36:47Z centler $
      subroutine newtonsub(j,iflag,ftol1,etol1,newton,ilud)

      include 'common_geo.inc'
      include 'common.inc'

       dimension funcs(ncomp),pd(ncomp,ncomp)
       dimension indexx(ncomp),ctemp(ncomp)
       dimension pdo(ncomp,ncomp),y(ncomp)

      integer newton
      integer ilud
      integer iflag

c***********************************************************************
c  START OF NEWTON ITERATIONS                                          *
c   NOTE: CHANGES HOW TO DEAL WITH FIXED CONCENTRATION BOUNDARY POINTS *
c     newton = max. number of iterations                               *
c     ne     = iteration number                                        *
c     funcs  = initially the vector of function residuals (Maple),     *
c              then replaced by calculated change in concentrations in *
c               each newton step                                       *
c     fxmax  = maximum residual in the current time step & position    *
c                 (residual x delt).                                   *
c     pd     = 2D array representing the Jacobian matrix (Maple)       *
c              evaluated for concentrations before transport           *
c     pdo    = LU decomposed version of pd                             *
c     indexx = vector of size ncomp to used by ludcmp and lubksb       *
c              (row permutations)                                      *
c     det    = parity of the number of row inerchanges                 *
c     y      = temporary vector of residuals/change in conc., see funcs*
c     ctemp  = temporary vector of concentrations                      *
c     relax  =                                                         *
c     cmplim =                                                         *
c     r2     =                                                         *
c     errmax = maximum calculated change in concentration for all      *
c              components in this position and time step               *
c     residual = routine containing function residuals (Maple)         *
c     jacobian = routine containing the entries of the Jacobian matrix *
c     ludcmp   = LU decomposition (see Numerical recipes, pp 39-51)    *
c     lubksb   = LU backsubstitution                                   *
c     mprove   = solution "improver" (see Numerical Recipes, p.51)     *
c     j      = spatial node                                            *
c     iflag  = warning flag, keeps track of neg. conc a.s.o            *
c     etol1  = acceptable absolute deviation from root                 *
c     ftol1  = acceptable relative deviations from root                *
c***********************************************************************

        iflag = 0

        do ne = 1,newton
          call residual(funcs,j)
          fxmax = 0.d0
          do i = 1,ncomp
            funcs(i) = -funcs(i)
            if(dabs(funcs(i)).gt.fxmax) then
              fxmax = dabs(funcs(i))*delt
            endif
          end do
          call jacobian(pd,j)


c     ****alternative for inverting the jacobian
c     igauss=1: gauss jordan, else: LUDecomposition
c     in the gauss jordan section, I also compute now the
c     eigenvalues and eigenvectors for the jacobian
c     this is for testing purposes and has nothing to do
c     with the simulations.
c     feel free to comment eigen.f out
c     GAUSS-JORDAN with full pivoting
c     in: jacobian, size of jacobian,physical size,-f,
c          nr.rhs-vectors,physical size)
c     out: inverse jacobian,size,size,solution (i.e. deltaC),1,1
      igauss=0 !1
      if(igauss.eq.1)then
!          call eigen(pd,ncomp)
           call gaussj(pd,ncomp,ncomp,funcs,1,1)
          do i = 1,ncomp
            ctemp(i) = sp(i,j)
          end do
      else
c     LUD
          do i = 1,ncomp
            do i2 = 1,ncomp
              pdo(i,i2) = pd(i,i2)
            end do
          end do

          do i = 1,ncomp
            y(i) = funcs(i)
          end do

c     LUD choice
       if(ilud.eq.1)then
           call ludcmp(pdo,ncomp,ncomp,indexx,det)
           call lubksb(pdo,ncomp,ncomp,indexx,y)
         else
c      If the linker complains about missing dgesv_ symbol, the next
c      five lines could be removed. But double check that the function
c      is never called ...
           call DGESV(ncomp,1,pdo,ncomp,indexx,y,ncomp,info)
           if(info.ne.0)then
                if(info.lt.0)write(*,*)'illegal value in LUD', info
                if(info.gt.0)write(*,*)'singular matrix in LUD'
           endif
       endif

          if (ne.gt.8) then
            call mprove(pd,pdo,ncomp,ncomp,indexx,funcs,y)
            call mprove(pd,pdo,ncomp,ncomp,indexx,funcs,y)
          end if

          do i = 1,ncomp
            funcs(i) = y(i)
            ctemp(i) = sp(i,j)
          end do
      endif

          relax= 1.0d0

          do i=1,ncomp
            r2 = 1.d0
            if(ctemp(i).gt.0.d0) then
              cmplim=0.1d0*ctemp(i)
              func_mag = dabs(funcs(i))
              if (func_mag.gt.cmplim) then
                r2= dabs(cmplim/funcs(i))
              else
                r2=1.d0
              end if
            endif
            sp(i,j) =  sp(i,j) + funcs(i)*r2

c***********************************************************************
c  if fixed concentration at boundary, DO NOT adjust but preserve      *
c  boundary condition also during the root finding                     *
c  avoid newton-looping due to fixed conc. boundary condition -> dC=0  *
c***********************************************************************
            if ((j.eq.1).and.(ibc(i,1).eq.0)) then
              sp(i,j) = spb(i,1) ! preserve boundary
             funcs(i) = 0.0d0
           end if
            if ((j.eq.nx).and.(ibc(i,2).eq.0)) then
              sp(i,j) = spb(i,2) ! preserve boundary
             funcs(i) = 0.0d0
           end if

            if (sp(i,j).lt.0.d0) iflag = 1

          end do

          errmax = 0.d0
          do i = 1,ncomp
            if(dabs(funcs(i)).gt.errmax) then
              errmax = dabs(funcs(i))
            endif
          end do
          if(errmax.lt.etol1 .or. fxmax.lt.ftol1) return

        end do  ! end of newton loop

        ! mark flag if exceeding newton interations
        if (iflag.eq.0) iflag=2
        if (iflag.eq.1) iflag=3

        return
        end
