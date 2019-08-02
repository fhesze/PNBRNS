c      
c     SUBROUTINE residual
c      
      subroutine residual(funcs,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension funcs(ncomp)
        call switches(j)
          funcs(1) = -1/radi*mumax*robac*sp(1,j)/(km+sp(1,j))/250-sp(1,j
     +)/delt+spold(1,j)/delt
          funcs(2) = -0.4934511124D1*Diff/radi**2*(sp(2,j)+km+0.40530864
     +15D-3*radi*mumax*robac/Diff)*(0.1D1-(0.1D1-0.1621234566D-2*sp(2,j)
     +*radi*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*radi*mumax*roba
     +c/Diff)**0.2D1)**0.5D0)-sp(2,j)/delt+spold(2,j)/delt
          funcs(3) = -1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))*kinhib/(ki
     +nhib+sp(3,j))/250-sp(3,j)/delt+spold(3,j)/delt
          funcs(4) = sp(4,j)-spold(4,j)
          funcs(5) = sp(5,j)-spold(5,j)
      end
