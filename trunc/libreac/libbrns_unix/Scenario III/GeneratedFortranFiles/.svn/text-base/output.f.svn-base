c      
c     SUBROUTINE out
c      
      subroutine out(j,nt,time,depth,v_out,v_int)
        include 'common_geo.inc'
        include 'common.inc'
        real*8 time
        if (nt.eq.1.and.j.eq.1) then
        v_out = 0.1D-1
        v_int = 0.2D-1
        endif
        if (time.le.v_out.and.v_out.lt.time+delt) then
        if (j.eq.nx) then
        v_out = v_out+v_int
        endif
        endif
        if (time.eq.endt) then
        endif
      end
