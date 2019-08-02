c      
c     SUBROUTINE jacobian
c      
      subroutine jacobian(pd,j)
        include 'common_geo.inc'
        include 'common.inc'
        dimension pd(ncomp,ncomp)
        call switches(j)
         pd(4,2) = 0
         pd(1,4) = 0
         pd(3,1) = 0
         pd(1,2) = 0
         pd(2,2) = -0.4934511124D1*Diff/radi**2*(0.1D1-(0.1D1-0.16212345
     +66D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*ra
     +di*mumax*robac/Diff)**0.2D1)**0.5D0)+0.2467255562D1*Diff/radi**2*(
     +sp(2,j)+km+0.4053086415D-3*radi*mumax*robac/Diff)/(0.1D1-0.1621234
     +566D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.4053086415D-3*r
     +adi*mumax*robac/Diff)**0.2D1)**0.5D0*(-0.1621234566D-2*radi*mumax*
     +robac/Diff/(sp(2,j)+km+0.4053086415D-3*radi*mumax*robac/Diff)**0.2
     +D1+0.3242469132D-2*sp(2,j)*radi*mumax*robac/Diff/(sp(2,j)+km+0.405
     +3086415D-3*radi*mumax*robac/Diff)**0.3D1)-1/delt
         pd(4,3) = 0
         pd(5,1) = 0
         pd(1,5) = 0
         pd(3,2) = 0
         pd(2,3) = 0
         pd(4,4) = 1
         pd(5,2) = 0
         pd(2,4) = 0
         pd(4,5) = 0
         pd(3,3) = -1/radi*mumax*robac/(km+sp(3,j))*kinhib/(kinhib+sp(3,
     +j))/250+1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))**2*kinhib/(kinhib+
     +sp(3,j))/250+1/radi*mumax*robac*sp(3,j)/(km+sp(3,j))*kinhib/(kinhi
     +b+sp(3,j))**2/250-1/delt
         pd(1,1) = -1/radi*mumax*robac/(km+sp(1,j))/250+1/radi*mumax*rob
     +ac*sp(1,j)/(km+sp(1,j))**2/250-1/delt
         pd(2,1) = 0
         pd(5,3) = 0
         pd(2,5) = 0
         pd(3,4) = 0
         pd(5,4) = 0
         pd(3,5) = 0
         pd(5,5) = 1
         pd(4,1) = 0
         pd(1,3) = 0
      end
