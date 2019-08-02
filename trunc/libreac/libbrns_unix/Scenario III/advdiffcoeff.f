c    $Id: advdiffcoeff.f 16 2007-10-18 12:36:47Z centler $
c     $Id: advdiffcoeff.f 16 2007-10-18 12:36:47Z centler $
      subroutine advdiffcoeff()

      include 'common_geo.inc'
      include 'common.inc'

c***********************************************************************
! define your advection velocities and diffusion coefficients.         *
! advection and dispersion coefficients are defined at faces and nodes,*
! i.e. 1,2,...,nx  (the values @ nodes are only needed for the flux or *
! gradient boundary conditions)                                        *
! w acts on solids and solutes                                         *
! q acts only on solutes                                               *
! Dmol acts only on solutes (can be overwritten)                       *
! Disp acts only on solutes                                            *
! Db acts on solids and solutes                                        *
! depth can be referred to as x(j)                                     *
!                                                                      *
! disp is the total diffusion coefficient for a species k at face j    *
! vd is the advection velocity of a dissolved species at face j        *
! vp is the advection velocity of a solid species at face j            *
!                                                                      *
! CM, July 2002                                                        *
c***********************************************************************

      real*8 dsol,dmix,db,fac_sal,vq
      integer isolid

       do j=1,nx
        vd(j) = 0.0d0
        vs(j) = 0.0d0
        if (iw.eq.0) then ! constant w
           w = w0
        else
!         user defined w
!         w = ...
        end if

        if (iq.eq.0) then ! constant q
          vq = q0/(por(j)*area(j))
        else
!         user defined v
!         vq = ...
        end if

         vd(j) = vd(j) + w + vq
         vs(j) = vs(j) + w

        dmix = aL*dabs(vq)          ! dispersion

        if(idb.eq.0) then               ! bioturbation
          db = db0
         else
           ! user defined bioturbation coefficient
!         db = ...
c          db = db0/2.0d0*erfcc((x(j)-20.0d0)/4.0d0)
           db = db0/2.0d0*erfcc((x(j)-25.0d0)/4.0d0)
c          db = db0/2.0d0*erfcc((x(j)-10.0d0)/4.0d0)
c         db=0.1d0
c         if (x(j).gt.15.0d0) db=0.d0
        end if

         do k = 1,ncomp
          disp(k,j) = 0.0d0

           fac_sal = (9.5d-01 - 1.d-03 * salin)
          dsol = fac_sal*dsol_0(k)*(1.d0+f_T(k)*t_celsius)! Dmol(T,S)
           dsol = dsol/(1.0d0-dlog(por(j)*por(j)))       ! in situ Dmol

           isolid=0
           call issolid(k,isolid)     ! set isolid to 1 if k is a solid
          if (isolid.eq.1) then
            disp(k,j) = db
            if (dsol.gt.0.0d0) then
              write(*,*) 'molecular diffusion for solid? ok...'
              disp(k,j) = disp(k,j)+dsol
              pause 'hit any key to continue'
            end if
          else
             disp(k,j) = dmix + dsol + db
           end if
         enddo

       enddo

c***********************************************************************
c   additional properties outside the domain                           *
c   required for flux boundary conditions, done by linear extrapolation*
c***********************************************************************
       do ii=1,2
        if (ii.eq.1) then
          i = 1
          id = -1
        end if
        if (ii.eq.2) then
          i = nx
          id = 1
        end if
         vs(i+id)=vs(i)+(vs(i-id)-vs(i))/(x(i-id)-x(i))*(x(i+id)-x(i))
         vd(i+id)=vd(i)+(vd(i-id)-vd(i))/(x(i-id)-x(i))*(x(i+id)-x(i))
        do k=1,ncomp
           disp(k,i+id) = disp(k,i) +
     +         (disp(k,i-id)-disp(k,i))/(x(i-id)-x(i))*(x(i+id)-x(i))
        end do
      end do

c***********************************************************************
c overwrite ghost points here below if desired                         *
c the depth of the ghost points is set in gridsetup.f                  *
c***********************************************************************
!      vs(0) =
!      vd(0) =
!      vs(nx+1) =
!      vd(nx+1) =
!      do k=1,ncomp
!        disp(k,0) =
!        disp(k,nx+1) =
!      end do


c***********************************************************************
c   check of physics: 0>=disp                                          *
c***********************************************************************
       do ii=0,nx+1
        do k=1,ncomp
           if (disp(k,ii).lt.0.0d0) then
            write(*,*) 'disp<0->set to 0', ii, x(ii)
            disp(k,ii)=0.0d0
             !stop
          end if
         end do
      end do

      end


c***********************************************************************
c   complimentary error function as function available                 *
c   based on Chebyshev fitting, taken from Press et al.                *
c   promoted to double precision                                       *
c***********************************************************************
       FUNCTION erfcc(x)
       REAL*8 erfcc,x
       REAL*8 t,z
       z=dabs(x)
       t=1./(1.d0+0.5d0*z)
       erfcc=t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     *(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     *(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
       if (x.lt.0.d0) erfcc=2.d0-erfcc
       return
       END


