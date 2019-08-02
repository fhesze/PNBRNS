c      
c     SUBROUTINE ssrates
c      
      subroutine ssrates(rat,drdc,isp,j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            if (isp.eq.1) then
            rat = -1/radi*mumax*robac*sp(1,j)/(km+sp(1,j))/250
            drdc = -1/radi*mumax*robac/(km+sp(1,j))/250+1/radi*mumax*rob
     +ac*sp(1,j)/(km+sp(1,j))**2/250
            endif
            if (isp.eq.2) then
            rat = -0.4934511124D1*Diff/radi**2*(sp(2,j)+km+0.4053086415D
     +-3*radi*mumax*robac/Diff)*(0.1D1-(0.1D1-0.1621234566D-2*sp(2,j)*ra
     +di*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*radi*mumax*robac/D
     +iff)**0.2D1)**0.5D0)
            drdc = -0.4934511124D1*Diff/radi**2*(0.1D1-(0.1D1-0.16212345
     +66D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*ra
     +di*mumax*robac/Diff)**0.2D1)**0.5D0)+0.2467255562D1*Diff/radi**2*(
     +sp(2,j)+km+0.4053086415D-3*radi*mumax*robac/Diff)/(0.1D1-0.1621234
     +566D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*r
     +adi*mumax*robac/Diff)**0.2D1)**0.5D0*(-0.1621234566D-2*radi*mumax*
     +robac/Diff/(sp(2,j)+km+0.4053086415D-3*radi*mumax*robac/Diff)**0.2
     +D1+0.3242469132D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.405
     +3086415D-3*radi*mumax*robac/Diff)**0.3D1)
            endif
            if (isp.eq.3) then
            rat = -1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))*kinhib/(kinhi
     +b+sp(3,j))/250
            drdc = -1/radi*mumax*robac/(km+sp(3,j))*kinhib/(kinhib+sp(3,
     +j))/250+1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))**2*kinhib/(kinhib+
     +sp(3,j))/250+1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))*kinhib/(kinhi
     +b+sp(3,j))**2/250
            endif
            if (isp.eq.4) then
            rat = 0.D0
            drdc = 0
            endif
            if (isp.eq.5) then
            rat = 0.D0
            drdc = 0
            endif
      end
