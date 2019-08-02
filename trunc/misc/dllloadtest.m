% Initialization

loadlibrary('brnsdll.dll', 'brnsdll.h')

% Dummy variables; they must be there, but have no effect
fixedConcentrationBoundary=[0 0] % must have entry for every species!
waterSaturation=0
porosity=0
pos_x=0
pos_y=0
pos_z=0
fcb_ptr=libpointer('int32Ptr', fixedConcentrationBoundary)
p_ptr=libpointer('doublePtr', porosity)
ws_ptr=libpointer('doublePtr', waterSaturation)
pos_x_ptr=libpointer('doublePtr', pos_x)
pos_y_ptr=libpointer('doublePtr', pos_y)
pos_z_ptr=libpointer('doublePtr', pos_z)

returnValue=-1
rv_ptr=libpointer('int32Ptr', returnValue)

numberOfSpecies=2
timeStep=0.1

timeStep_ptr=libpointer('doublePtr', timeStep)
numberOfSpecies_ptr=libpointer('int32Ptr', numberOfSpecies)

parameterVector=[]
parameterVector_ptr=libpointer('doublePtr', parameterVector)

% Reactive Step in time loop

ConcBeforeTransport=[0.5 0.5]
ConcAfterTransport =[0.6 0.6]

% Coupling Starts

cat_ptr=libpointer('doublePtr', ConcAfterTransport)
cbt_ptr=libpointer('doublePtr', ConcBeforeTransport) % used as initial guess for newton iteration

% Calling BRNS

calllib('brnsdll','invokebrns',cat_ptr, cbt_ptr, cat_ptr, numberOfSpecis_ptr, timeStep_ptr, fcb_ptr, rv_ptr, pos_x_ptr, pos_y_ptr, pos_z_ptr, p_ptr, ws_ptr, parameterVector_ptr)

% Retrieving Data

returnValue=get(rv_ptr, 'Value')
ConcBeforeTransport=get(cat_ptr, 'Value') 

% CHECK returnValue:
% 0: Computation ok
% 1: Negative concentrations occured
% 2: Exceeding maximum newton iteration
% 3: 1 and 2

% Coupling Ends

% Before program exit

unloadlibrary('brnsdll')
