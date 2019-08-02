c    $Id: porarea.f 16 2007-10-18 12:36:47Z centler $
      subroutine porarea()

        include 'common_geo.inc'
        include 'common.inc'

c***********************************************************************
! define porosity and cross sectional area                             *
! they are both defined at nodes and faces, i.e. j=1,2,..nx-1,nx       *
! depth can be referred to as x(j)                                     *
! note the ghost point below; see also comments in gridsetup.f         *
!                                                                      *
! CM, July 2002                                                        *
c***********************************************************************
        real*8 dist(nx)

      do j=1,nx
        if (ipor.eq.0) then     ! constant porosity
           por(j) = por0
        else
!         user defined
!         por(j) =
        endif

        if (iarea.eq.0) then     ! constant cross section area
           area(j) = area0
        else
!         user defined
!         area(j) =
           if (j.eq.1) then
              dist(1) = 0.0d0
          else
             dist(j) = dist(j-1) + (x(j+1)-x(j-1))/2000.0d0
          end if
           area(j) = 151350.0d0 * 10.0d0**(-0.02d0*dist(j)) ! 30000.0d0
        endif
       end do

c***********************************************************************
c   additional por and area definition outside the domain              *
c   required for flux boundary conditions                              *
c***********************************************************************
       por(0)     = por(1)   - (por(2)-por(1))
       area(0)    = area(1)  - (area(2)-area(1))
       por(nx+1)  = por(nx)  + (por(nx)-por(nx-1))
       area(nx+1) = area(nx) + (area(nx)-area(nx-1))


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
       por(i+id)=por(i)+(por(i-id)-por(i))/(x(i-id)-x(i))*(x(i+id)-x(i))
         area(i+id)=area(i)+
     +                (area(i-id)-area(i))/(x(i-id)-x(i))*(x(i+id)-x(i))
       end do

c***********************************************************************
c overwrite ghost points here below if desired                         *
c the depth of the ghost points is set in gridsetup.f                  *
c***********************************************************************
!      por(0) =
!      area(0) =
!      por(nx+1) =
!      area(nx+1) =

c***********************************************************************
c   check of physics: 0>=por>=1, area >0                               *
c***********************************************************************
       do ii=0,nx+1
        if (por(ii).gt.1.0d0) then
          write(*,*) 'por>1:', ii, x(ii)
          stop
        end if
        if (por(ii).le.0.0d0) then
          write(*,*) 'por<=0:', ii, x(ii)
          stop
        end if
         if (area(ii).le.0.0d0) then
          write(*,*) 'area<=0:', ii, x(ii)
          stop
        end if
      end do

      end
