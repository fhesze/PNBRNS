c    $Id: gridsetup.f 16 2007-10-18 12:36:47Z centler $
      subroutine gridsetup()

c***********************************************************************
! define your grid                                                     *
! i.e. the depths of nodes and faces (1,2...nx-1,nx)                   *
! - note that the nodes are recommened to be in the middle of the box  *
! (because the concentration at the interface is done by linear        *
! interpolation, the concentration in at node i is in my opinion an    *
! average for the box. Thus, linear interpolation to the faces should  *
! originate from the middle of the box. Therefore, the node should be  *
! in the middle of two faces. If you agree, this can be done by        *
! defining only the depth of the nodes, i.e. fix node 1 at a depth of 0*
! and then step from j=3 in steps of 2 (i.e. only define nodes).       *
! Then calculate the faces as x(j-1) = x(j)-x(j-2) )                   *
! - If you implement a user-defined grid, you might also might want to *
! check that x(nx)=depthmax                                            *
! - It is also recommended to only make the grid vary somewhat slowly  *
! because formally you loose accuracy. For suggestions see Noye's book *
! - If you make funky grids, the best thing to do is to make D, v, por *
! and area, which must be provided in the file advdiffcoeff.f and      *
! porarea.f, a function of depth (distance) and refer to depth as a    *
! function of x(j). To check the input, the grid and properties are    *
! echoed in echogrid.dat                                               *
!                                                                      *
! for flux or gradient boundary conditions, ghost points outside the   *
! domain are required. This is done automatically (see bottom).        *
! However, if you know the properties (D, por a.s.o) at a specific     *
! point outside the domain, you might want to specify these position   *
! explicitly; the ghost points have indices 0 and nx+1.                *
!                                                                      *
! CM, July 2002                                                        *
c***********************************************************************

       include 'common_geo.inc'
        include 'common.inc'

       x(1) = 0.0d0
      if (igrid.eq.0) then ! even grid
         dx = depthmax/dfloat(nx-1)
        do j=2, nx
          x(j) = x(j-1) + dx
        end do
      else


        do j=3,nx,2

c          if (x(j-2).lt.29.99d0) then
c            x(j) = x(j-2) + 0.3d0
c          else if (x(j-2).lt.47.99d0) then
c            x(j) = x(j-2) + 0.6d0
c          else if (x(j-2).lt.53.99d0) then
c           x(j) = x(j-2) + 1.2d0
c         else if (x(j-2).lt.65.99d0) then
c           x(j) = x(j-2) + 2.4d0
c         else if (x(j-2).lt.89.99d0) then
c            x(j) = x(j-2) + 4.8d0
c         else
c            x(j) = x(j-2) + 9.d0
c          end if
           if (x(j-2).lt.3.0d0) then
                x(j)=x(j-2)+0.1d0
           else if (x(j-2).lt.5.0d0) then
                x(j)=x(j-2)+0.2d0
           else if (x(j-2).lt.10.0d0) then
                x(j)=x(j-2)+0.25d0
           else if (x(j-2).lt.20.0d0) then
                x(j)=x(j-2)+0.5d0
           else if (x(j-2).lt.50.0d0) then
                x(j)=x(j-2)+1.0d0
           else
                x(j)=x(j-2)+2.5d0
           end if

            x(j-1) = x(j-2) + (x(j)-x(j-2))/2.d0
        end do

        if (x(nx).ne.depthmax) then
          write(*,*) 'last node = depthmax ?', x(nx), depthmax
        end if
      end if


c***********************************************************************
c   additional depth definition outside the domain                     *
c   required for flux boundary conditions                              *
c   if you know D, por and area at locations outside the domain, you   *
c   might want to set the correct depths here, rather than use a linear*
c   extrapolation. make sure that x(0)<x(1), x(nx+1)>x(nx)             *
c***********************************************************************
       x(0)  = x(1) - (x(2)-x(1))
       x(nx+1) = x(nx)   + (x(nx)-x(nx-1))

      ! x(0) =
      ! x(nx+1) =

      end
