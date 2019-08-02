# _________________________________________________________________________

# AUTOMATIC CODE GENERATOR (ACG) FOR CONSTRUCTING USER-DEFINED
# BIOGEOCHEMICAL REACTION NETWORKS
# 
# version 1.2
# COPYRIGHT (c) 2001 P.A.G. Regnier 
# All Rights Reserved
# 
# Research Unit on Biogeochemical Systems Dynamics
# Department of Geochemistry, Utrecht University, 
# The Netherlands
# __________________________________________________________________________
# __________________________________________________________________________
# INPUT TYPES
# 
# OOO : Sections that should be modified by the user
# OOO : Sections that should NOT be modified by the user

# OOO : comments 
# OOO : Maple input
# OOO : Maple output (appears only after you have executed the spreadsheet) 
# OOO  : Maple input entries that have to be specified by the user 
# 
# WWW  : Hyperlink to the Knowledge Book
# 
# _________________________________________________________________________ 

#  Maple specific info
restart ;
#with(Spread) :
 precision := double :
#  Summary and governing equations
#  Caveats 
# - for real numbers you should add a point after the number.
# - make sure all units match. 
# e.g. rates and dXdt, or flux boundary conditions and concentrations
# for the conversion of solid to solute units, 
# one may define temporary variables that can be used in the rate laws below:
# s_dens := 2.5; # solid density in [g/cm_solid^3]
# sd := 1000. * s_dens * (1. - por(j)) / por(j); 
# the factor 1.d03 converts cm^3 to liter, e.g. [g/cm^3] to [g/l]. 
# note that you need to refer to porosity exactly as por(j).
#  Reaction Network - Size and Variables
#  Size of reaction network
# nsolids : number of solid species 
# ndissolved : number of dissolved species 
# ncompo : total number of species 
# nreactions : total number of reactions (including equilibrium rxns) 
# neqrxns : number of equilibrium reaction
ncompo := 4 ;
nreactions := 3 ;
neqrxns := 0 ;
#  List of variables
# variables: list of variables to model. example: 
# listsolids: species number which is a SOLID species. 
# note: all other variables are temporary and are NOT parsed to the ACG
# 
# Example: 
# variables:=[O2, so4, MnOx, FeOx, hco3, co3, hplus, hs];
# listsolids: = [3,4];
 variables := [cbioeff,cbio,corg,tracer] ;
#  Biogeochemistry - Rate laws
# Definition of kinetic rate laws
# rate.i : array of rates. 
# - For equilibrium rate expression, a kinetic rate MUST be specified as well. It will be overwritten in the equilibrium section below, but you need it as space holder and stoichiometry. Furthermore, the steadystate module uses detailed balancing method with fast kinetics. Therefore, in the example below, kf will have to be defined as a rather large number and kb = kf*Keq
# note: all other variables are temporary and are NOT parsed to the ACG
# - conditional statements: if a rate law depends on a conditional statement you need to make use of the subroutine switches.f. Example: dissolution (Rd) is only to take place at undersaturation, thus  Rd= f(saturation). If saturation > 1, Rd>0, else Rd=0. This canbe implemented as Rd:=k*H1*saturation, where H1 is toggled between 0 and 1. Rather than giving the condition here in maple, for now you need to do this in "switches.f", where you program the conditions for H1, e.g. H1=0, If (A*B/K>1) then H1 = 1.
#  
# example: 
# rate1 := 1000.*O2*hs; # rate law for 2O2 + HS -> SO4 + Hplus
# rate2 := kf*hplus*co3 - kb*hco3; # kinetic rate law for HCO3 = CO3 + Hplus (equilibrium)
# 
#  Primary redox reactions WWW

#  Sulfide Re-oxidation reactions

#k Best,My interpretation, Dhyd and Av from Cylinder
#Av:=2/radi;
#Dhyd:=2*radi;
#kexch:=(3.1415**2.0)*Diff*Av/(4.0*Dhyd);

#k Best,Martin interpretation, rhyd=radi, Av from experiment=2 times of 160mic
Av:=4/radi;
kexch:=(3.1415**2.0)*Diff*Av/(4.0*radi);

#k Best
#kexch:=(3.1415**2.0)*Diff/(4.0*(radi**2.0));

#robac:=massbac/surfexp;
#increase the reactivity of surface by factor of 2
#robac:=2*massbac/surfexp;

#qmax:=2.0*mumax*robac/radi/1000;
qmax:=Av*mumax*robac/1000;


term1:=qmax/kexch;
term21:=cbioeff+km+term1;
term22:=cbio+km+term1;

rate1:=kexch/2.0*term21*(1.0-(1.0-4.0*cbioeff*term1/(term21**2.0))**0.5);
rate2:=kexch/2.0*term22*(1.0-(1.0-4.0*cbio*term1/(term22**2.0))**0.5);

rate3 := qmax*corg/(km+corg) ;

#Kfac:=1+corg/km+4*mumax*robac*radi/(3.14)^2/Diff/km;
#rate1:= (3.14)^2*Diff*km/radi^2*Kfac-(1-sqrt(1-8*corg/km*mumax*robac*radi/(3.14)^2/Diff/km/Kfac^2)) ;
#rate1:=radi*corg;
# rate2:= muo*baco*c_bioav/(kmoc+c_bioav)*o2/(kmo+o2) ; #*(1-(baco+bacn+bacs)/baccmax);
# rate3:= deco*baco*(1-bacmin/baco) ;
# rate4:= mus*bacs*c_bioav/(kmsc+c_bioav)*so4/(kms+so4)*kmoinh/(kmoinh+o2); #*(1-(baco+bacn+bacs)/baccmax); 
# rate5:= decs*bacs*(1-bacmin/bacs) ;

#  Biogeochemistry - Stoichiometry

# Stoichiometry of the biogeochemical reactions
# d.sp.dt : rates of change of sp due to the sum of biogeochemical reactions
# note that rateX must be referred to as rX  
#  
# example:
# dO2dt := -2*r1;
# dhco3dt = -r2; 
s_dens := 2.5 ;
#SD := 1.0e3 * s_dens * (1.0 - por(j)) / por(j);
# 

dcbioeffdt:=-r1;
dcbiodt:= -r2;
dcorgdt:=-r3;
dtracerdt:=0.0;



#  Biogeochemistry - Equilibria
# Specification of equilibrium constraints
# eqrxnsId : set of kinetic reactions which are overuled by a thermodynamic constraint 
# equilibriumseqns[i] : Equilibrium constraint for reaction i  
# 
# example:
# eqrxnID := [r2,rX];
# equilibriumeqns[1] := hplus*co3 - Keq*hco3;
# equilibriumeqns[2] := ...;
eqrxnId := [] ;
#equilibriumeqns[1] := b - k3*c ;
#  Biogeochemistry - Parameters
# Values of rates constants and parameters 
# In this section, all parameters defined in section 'Rate laws' should be defined.
# nparam: number of parameters to define 
# The list is given by bio_name; the values collected in bio_val.
# note that for double precision, 10 should be written as 10.
# 
# example:
# nparam:=4;
# bio_name:=[kmo2hs,kf,kb,Keq];
# vkf :=1.0*10^(5);
# vKeq:=1.0*10^(-10.4);
# vkb :=vkf*vKeq;
# bio_val:=[1000.,vkf,vkb,vKeq];
# 
# 
nparam := 5 ;
#bio_name := [km,mumax,robac,Diff,radi] ;
#bio_val :=[231.0e-9,3.3e-8,4.1,6.0e-10,0.0] ; 
bio_name := [km,mumax,Diff,radi,robac] ;
bio_val :=[231.0e-9,0.326e-9,6.0e-10,0.0,0.0] ; 



# Parameters
# Parameters can be used in the rate equations for parameters whose values are passed to the reactive solver each time it is called. Specify in nparameters, how many parameters you need and in parameterslist the names of the parameters. The parameter names must also appear in bio_name and must be assigned a dummy value there. When calling the library, use the order as listed here to pass the values. It is your responsibility to supply the correct number of values!
nparameters := 2 ;
parameterslist := [radi,robac] ;
# 
#   
# Switches
# Switches can be used in the rate equations. Specify in nswitches, how many switches are in use, name them and define the switch expressions. The switch names must also appear in bio_name and must be assigned a dummy value there. The switch equals 1 if the switch expression is >0, 0 otherwise. To reference the coordinates in the domain, use x_pos, y_pos and z_pos.
nswitches := 0 ;
switchlist := [sw1,sw2,sw3] ;
switchcrit := [(bio*dissb-0.25),dissc,dissc] ;
# 
#  Maple specific info
# dir_f: directory where the FORTRAN routines and Maple spread.m files are parsed
# format Mac: "Macinthosh HD:UU:...:code"
# format PC: "C:\\maple\\...\\code"
# WAS: dir_f := "C:\\Dokumente und Einstellungen\\centler\\Desktop\\Labor\\Simulations": 
# currentdir(dir_f):
# save "spread.m" ;
dir_f := "C:\\Program Files\\BRNSPackage" :
currentdir(dir_f) :
parse(sprintf("save %q,\"spread.m\";",anames()), statement) ;

"now execute processor - make sure the directories are set correctly";
# 
# ACG
