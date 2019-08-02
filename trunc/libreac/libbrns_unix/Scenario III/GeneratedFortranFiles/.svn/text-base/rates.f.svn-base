c      
c     SUBROUTINE rates
c      
      subroutine rates(j)
        include 'common_geo.inc'
        include 'common.inc'
        call switches(j)
            r(1,j) = 1/radi*mumax*robac*sp(1,j)/(km+sp(1,j))/250
            r(2,j) = 0.4934511E1*Diff/radi**2*(sp(2,j)+km+0.4053086E-3*r
     +adi*mumax*robac/Diff)*(0.1E1-(0.1E1-0.1621235E-2*sp(2,j)*radi*muma
     +x*robac/Diff/(sp(2,j)+km+0.4053086E-3*radi*mumax*robac/Diff)**0.2E
     +1)**0.5E0)
            r(3,j) = 1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))*kinhib/(kin
     +hib+sp(3,j))/250
            r(4,j) = 0.9869022E1*Diff/radi**2*(sp(4,j)-sp(5,j))
            r(5,j) = 1/radi*mumax*robac*sp(5,j)/(sp(5,j)+km)*kinhib/(kin
     +hib+sp(5,j))/250
      end
