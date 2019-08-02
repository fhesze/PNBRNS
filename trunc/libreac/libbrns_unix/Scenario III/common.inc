c    $Id: common.inc 16 2007-10-18 12:36:47Z centler $
         common/tsp/sp,co,spold,disp,dsol_0,f_T,r
!          common/tsp2/por,area,x,vd,vs,spb,spi
           common/tsp2/por,area,x,vd,vs,spb,spi,gm,ac
           common/coordinates/x_pos,y_pos,z_pos
         dimension sp(ncomp,nx),co(0:nx+1),spold(ncomp,nx),r(nreac,nx) 
           dimension disp(ncomp,0:nx+1),dsol_0(ncomp),f_T(ncomp)
           dimension por(0:nx+1),area(0:nx+1),x(0:nx+1),spb(ncomp,2)
           dimension vd(0:nx+1),vs(0:nx+1)
           common/boun/nssnow,ibc
           dimension ibc(ncomp,2)

           dimension gm((nx+1)*ncomp/2,(nx+1)*ncomp/2)
           dimension ac(ncomp*(nx+1)/2)
