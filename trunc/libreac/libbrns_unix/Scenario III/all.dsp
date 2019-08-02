# Microsoft Developer Studio Project File - Name="all" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=all - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "all.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "all.mak" CFG="all - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "all - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "all - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "all - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "all - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /cxml
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "all - Win32 Release"
# Name "all - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\advdiffcoeff.f
DEP_F90_ADVDI=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\basic.f
DEP_F90_BASIC=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\biogeo.f
DEP_F90_BIOGE=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\boundaries.f
DEP_F90_BOUND=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\diagenesis.f
DEP_F90_DIAGE=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\drivervalues.f
DEP_F90_DRIVE=\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\eigen.f
# End Source File
# Begin Source File

SOURCE=.\frngdb.f
# End Source File
# Begin Source File

SOURCE=.\funk.f
DEP_F90_FUNK_=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\GAMMP.F
# End Source File
# Begin Source File

SOURCE=.\gaussj.f
# End Source File
# Begin Source File

SOURCE=.\getdat.f
DEP_F90_GETDA=\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\getdelt.f
DEP_F90_GETDE=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\timestep.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\gridsetup.f
DEP_F90_GRIDS=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\initialcond.f
DEP_F90_INITI=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\interpolate.f
DEP_F90_INTER=\
	".\common.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\issolid.f
DEP_F90_ISSOL=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\jacobian.f
DEP_F90_JACOB=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\limits.f
DEP_F90_LIMIT=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\LUBKSB.F
# End Source File
# Begin Source File

SOURCE=.\LUDCMP.F
# End Source File
# Begin Source File

SOURCE=.\main.f
DEP_F90_MAIN_=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\molecular.f
DEP_F90_MOLEC=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\MPROVE.F
# End Source File
# Begin Source File

SOURCE=.\NEWT.F
DEP_F90_NEWT_=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\newtonsub.f
DEP_F90_NEWTO=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\notransport.f
DEP_F90_NOTRA=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\objf.f
DEP_F90_OBJF_=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\optima.f
DEP_F90_OPTIM=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\optimde.f
DEP_F90_OPTIMD=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\optimlm.f
DEP_F90_OPTIML=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\optimsa.f
DEP_F90_OPTIMS=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\output.f
DEP_F90_OUTPU=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\porarea.f
DEP_F90_PORAR=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\printdepth.f
DEP_F90_PRINT=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\rates.f
DEP_F90_RATES=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\residual.f
DEP_F90_RESID=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\ssrates.f
DEP_F90_SSRAT=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\steadystate.f
DEP_F90_STEAD=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\storedat.f
DEP_F90_STORE=\
	".\common.inc"\
	".\common_geo.inc"\
	".\common_meas.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\switches.f
DEP_F90_SWITC=\
	".\common.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\timestep.f
DEP_F90_TIMES=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	".\timestep.inc"\
	
# End Source File
# Begin Source File

SOURCE=".\transcoeff-MT.f"
DEP_F90_TRANS=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\transcoeff.f
DEP_F90_TRANSC=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\transferback.f
DEP_F90_TRANSF=\
	".\common_geo.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\transferfw.f
DEP_F90_TRANSFE=\
	".\common_geo.inc"\
	".\common_opt.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\transport.f
DEP_F90_TRANSP=\
	".\common.inc"\
	".\common_drive.inc"\
	".\common_geo.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\TRIDAG.F
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\blas_win32.lib
# End Source File
# Begin Source File

SOURCE=.\lapack_win32.lib
# End Source File
# End Target
# End Project
