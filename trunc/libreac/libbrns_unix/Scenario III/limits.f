c    $Id: limits.f 16 2007-10-18 12:36:47Z centler $
      subroutine limits(i,j)

       include 'common_geo.inc'
       include 'common.inc'

      integer j,i
            if ((sp(i,j).lt.1.0d30).and.(sp(i,j).gt.-1.0d30)) then
            ! do nothing. but if NAN enter in ...else...
             else
               write(*,*) 'out of bounds...reset to +/-1e30', i,j,
     &               sp(i,j)
#ifdef DLLWRITELOG
       write(25,*) 'out of bounds...reset to +/-1e30',i,j,sp(i,j)
#endif

               if(sp(i,j).gt.0.d0) then
                 sp(i,j) = 1.0d30
               else
                 if (sp(i,j).lt.0.d0) then
                   sp(i,j) = -1.0d30
                 else
                   sp(i,j) = 0.0d0
                 end if
               end if
c               write(*,*) j,i
             end if
      if (dabs(sp(i,j)).lt.1.0d-20) then
           if (sp(i,j).lt.0.0) then
                sp(i,j)=-1.0d-20
           else
                sp(i,j)=1.0d-20
           endif
      end if

      end

