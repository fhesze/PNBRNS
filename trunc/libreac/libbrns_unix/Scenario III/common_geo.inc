        implicit real*8(a-h,o-z) 
        parameter (nsolid=0,ndiss=5,ncomp=5,nreac=5,nx=5)
        common/kinetics/km,mumax,Diff,radi,robac,kinhib
        real*8 km,mumax,Diff,radi,robac,kinhib
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,endt
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,end
     +t
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic
        integer iq,iw,iDb,ipor,igrid,iarea,ic
