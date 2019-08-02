# $Id: makefile 16 2007-10-18 12:36:47Z centler $
#
#flag1=-O3 -qarch=pwr3 -qtune=pwr3 -qhot -p
#flag2=-c -O3 -qarch=pwr3 -qtune=pwr3 -qhot -p
#F77 = f90
CC = g++
F77 = gfortran
#F77 = g77
flag1= -O4 # -g #-ff90 -qzerosize
flag2=-c -O4 -x f77-cpp-input -I. -fPIC # -g #-ff90 -qzerosize  # FPIC PROBABLY DOESN'T HURT
#flag3= -O4 -fno-f2c -shared -s # -g #-ff90 -qzerosize		# FNO-F2C PROBABLY NOT REQUIRED
flag3= -L. -O4 -shared -s # -g #-ff90 -qzerosize
#flag4=-c -O4 -x f77-cpp-input -I. -fPIC # -g #-ff90 -qzerosize
flag4=-DRETURNRATES -c -O4 -x f77-cpp-input -I. -fPIC # -g #-ff90 -qzerosize 
# keep flag3 =-g
#flag3=-c -x f77-cpp-input -I. # -g #-ff90 -qzerosize
#flag4=-c -x f77-cpp-input -I. # -g #-ff90 -qzerosize
#flag1=-O4 -g
#flag2=-c -O4 -g
brns: advdiffcoeff.o basic.o biogeo.o boundaries.o diagenesis.o \
	drivervalues.o eigen.o frngdb.o funk.o GAMMP.o gaussj.o getdat.o getdelt.o \
	gridsetup.o initialcond.o interpolate.o issolid.o jacobian.o limits.o \
	LUBKSB.o LUDCMP.o main.o molecular.o MPROVE.o NEWT.o newtonsub.o\
	notransport.o objf.o optima.o optimde.o optimlm.o optimsa.o \
	output.o porarea.o printdepth.o rates.o residual.o \
	ssrates.o steadystate.o storedat.o switches.o \
	timestep.o transcoeff-MT.o transcoeff.o transferback.o transferfw.o transport.o TRIDAG.o printsvnversion_nosvn.o
	$(F77) $(flag1) advdiffcoeff.o basic.o biogeo.o boundaries.o diagenesis.o \
		drivervalues.o eigen.o frngdb.o funk.o GAMMP.o gaussj.o getdat.o getdelt.o \
		gridsetup.o initialcond.o interpolate.o issolid.o jacobian.o limits.o \
		LUBKSB.o LUDCMP.o main.o molecular.o MPROVE.o NEWT.o newtonsub.o\
		notransport.o objf.o optima.o optimde.o optimlm.o optimsa.o \
		output.o porarea.o printdepth.o rates.o residual.o \
		ssrates.o steadystate.o storedat.o switches.o \
		timestep.o transcoeff-MT.o transcoeff.o transferback.o transferfw.o transport.o TRIDAG.o printsvnversion_nosvn.o \
		-o brns -llapack -lblas

brns.so: basic.o biogeo.o boundaries.o drivervalues.o gaussj.o invokebrns.o \
	jacobian.o LUBKSB.o LUDCMP.o MPROVE.o newtonsub.o residual.o \
	switches.o parameters.o varporosity.o limits.o rates.o
	$(F77) $(flag3)  basic.o biogeo.o boundaries.o drivervalues.o \
	gaussj.o invokebrns.o jacobian.o LUBKSB.o LUDCMP.o MPROVE.o \
	newtonsub.o residual.o switches.o parameters.o varporosity.o \
	limits.o rates.o -o brns.so -llapack -lblas

CallingBrns: brns.so
	$(CC) -DUNIX_GNU -o CallingBrns CallingBRNS/CallingBRNS.cpp -ldl

advdiffcoeff.o: advdiffcoeff.f
	$(F77) $(flag2) advdiffcoeff.f
basic.o:basic.f 
	$(F77) $(flag4) basic.f
biogeo.o:biogeo.f
	$(F77) $(flag4) biogeo.f
boundaries.o:boundaries.f
	$(F77) $(flag4) boundaries.f
diagenesis.o:diagenesis.f
	$(F77) $(flag2) diagenesis.f
drivervalues.o:drivervalues.f
	$(F77) $(flag4) drivervalues.f
eigen.o:eigen.f
	$(F77) $(flag2) eigen.f
frngdb.o:frngdb.f
	$(F77) $(flag2) frngdb.f
funk.o:funk.f
	$(F77) $(flag2) funk.f
GAMMP.o:GAMMP.F
	$(F77) $(flag2) GAMMP.F
gaussj.o:gaussj.f
	$(F77) $(flag4) gaussj.f
getdat.o:getdat.f
	$(F77) $(flag2) getdat.f
getdelt.o:getdelt.f 
	$(F77) $(flag2) getdelt.f
gridsetup.o:gridsetup.f
	$(F77) $(flag2) gridsetup.f
initialcond.o:initialcond.f
	$(F77) $(flag2) initialcond.f
interpolate.o:interpolate.f
	$(F77) $(flag2) interpolate.f
issolid.o:issolid.f
	$(F77) $(flag2) issolid.f
jacobian.o:jacobian.f
	$(F77) $(flag4) jacobian.f
limits.o:limits.f
	$(F77) $(flag2) limits.f
LUBKSB.o:LUBKSB.F
	$(F77) $(flag4) LUBKSB.F
LUDCMP.o:LUDCMP.F
	$(F77) $(flag4) LUDCMP.F
main.o: main.f 
	$(F77) $(flag2) main.f
molecular.o:molecular.f
	$(F77) $(flag2) molecular.f
MPROVE.o:MPROVE.F
	$(F77) $(flag4) MPROVE.F
NEWT.o:NEWT.F
	$(F77) $(flag2) NEWT.F -I`pwd`
newtonsub.o:newtonsub.f
	$(F77) $(flag4) newtonsub.f
notransport.o:notransport.f
	$(F77) $(flag2) notransport.f
objf.o:objf.f
	$(F77) $(flag2) objf.f
optima.o:optima.f
	$(F77) $(flag2) optima.f
optimde.o:optimde.f 
	$(F77) $(flag2) optimde.f
optimlm.o:optimlm.f
	$(F77) $(flag2) optimlm.f
optimsa.o:optimsa.f
	$(F77) $(flag2) optimsa.f
output.o:output.f
	$(F77) $(flag2) output.f
porarea.o:porarea.f
	$(F77) $(flag2) porarea.f
printdepth.o:printdepth.f
	$(F77) $(flag2) printdepth.f
rates.o:rates.f
	$(F77) $(flag2) rates.f
residual.o:residual.f
	$(F77) $(flag4) residual.f
ssrates.o:ssrates.f
	$(F77) $(flag2) ssrates.f
steadystate.o:steadystate.f
	$(F77) $(flag2) steadystate.f
storedat.o:storedat.f
	$(F77) $(flag2) storedat.f
switches.o:switches.f
	$(F77) $(flag4) switches.f
parameters.o:BrnsDll/parameters.f
	( cd BrnsDll ; $(F77) $(flag4) -I.. parameters.f ; mv parameters.o .. )
varporosity.o:BrnsDll/varporosity.f
	( cd BrnsDll ; $(F77) $(flag4) -I.. varporosity.f ; mv varporosity.o .. )
timestep.o:timestep.f
	$(F77) $(flag2) timestep.f
transcoeff-MT.o:transcoeff-MT.f
	$(F77) $(flag2) transcoeff-MT.f
transcoeff.o:transcoeff.f
	$(F77) $(flag2) transcoeff.f
transferback.o:transferback.f
	$(F77) $(flag2) transferback.f
transferfw.o:transferfw.f
	$(F77) $(flag2) transferfw.f
transport.o:transport.f
	$(F77) $(flag2) transport.f
TRIDAG.o:TRIDAG.F
	$(F77) $(flag2) TRIDAG.F
printsvnversion_nosvn.o:printsvnversion_nosvn.f
	$(F77) $(flag2) printsvnversion_nosvn.f
invokebrns.o:BrnsDll/invokebrns.f
	( cd BrnsDll ; $(F77) $(flag4) -I.. invokebrns.f ; mv invokebrns.o .. )

clean:
	rm -f *.o
