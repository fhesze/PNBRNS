function spec = getReaction(FlowData, GeometryData, TransportCoeffs, spec)

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n           = GeometryData.GeometryCoeffs.xIncr;
m           = GeometryData.GeometryCoeffs.yIncr;
pore_len    = GeometryData.GeometryCoeffs.LengthOfPore;
node_no     = GeometryData.NodeData.NumberOfNodes;
nodes       = GeometryData.NodeData.NodesOfPores;
pore_pos    = GeometryData.PoreData.PorePosition;
poreXY      = GeometryData.PoreData.Pores;

% flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores 
pore_CS     = pi.*radi.^2;                                                 % pore cross section
pore_vol    = pi.*radi.^2.*pore_len; 

%%-- Transport Parameters -------------------------------------------------

c_spec_no   = 1;                                                           % number of chemical species 
b_spec_no   = 1;                                                           % number of bacterial species
spec_no     = c_spec_no + b_spec_no;                                       % total number of species

%%-- Time Parameters ------------------------------------------------------

t_end       = TransportCoeffs.t_end;
t_begin     = TransportCoeffs.t_0;
dt          = TransportCoeffs.dt;
% time_i      = 0;
% time_save   = 20; 

%%-------------------------------------------------------------------------

cNew        = spec.chem;
bac         = spec.bac;

%%--  BRNS file fixed params in lib-file ----------------------------------% this Av should be similar and somehow equal for hetrogenous networks of normally distributed values of radi

% mio_water   = 8.94e-4;                                                     % water viscosity
% Cond        = pi.*radi.^4./(8*mio_water*pore_len);                         % conductivity profile based on Haggen-P equation

%%-------------------------------------------------------------------------

pos_x           = 0;
pos_y           = 0;
pos_z           = 0;
posx_ptr        = libpointer('doublePtr', pos_x);
posy_ptr        = libpointer('doublePtr', pos_y);
posz_ptr        = libpointer('doublePtr', pos_z);

fixedConcBound  = zeros(spec_no);                                          % must have entry for every species!
waterSaturation = 0;
porosity        = 0;
fcb_ptr         = libpointer('int32Ptr', fixedConcBound);
p_ptr           = libpointer('doublePtr', porosity);
ws_ptr          = libpointer('doublePtr', waterSaturation);

returnValue     = 0;
rv_ptr          = libpointer('int32Ptr', returnValue);
nos_ptr         = libpointer('int32Ptr', spec_no);

timeStep        = dt;                                                      % for now we take same timestep from advection for reaction
ts_ptr          = libpointer('doublePtr', timeStep);

%%-------------------------------------------------------------------------

rand_num        = rand(pores_no,1);

%%-------------------------------------------------------------------------

path            = getLibPath('unix');
lib_file        = [path 'brns.so'];
lib_header      = [path 'brnsdll_param_new.h'];

[notfnd, warn]  = loadlibrary(lib_file, lib_header);

%%-------------------------------------------------------------------------

for pore_i=1:pores_no
        
    ConcAfterTransport  = [bac(pore_i) cNew(pore_i)];
    ConcBeforeTransport = ConcAfterTransport;                              % used as initial guess for newton iteration
    cOut(1:spec_no)     = 0.0;

    cat_ptr = libpointer('doublePtr', ConcAfterTransport);
    cbt_ptr = libpointer('doublePtr', ConcBeforeTransport);
    oc_ptr  = libpointer('doublePtr', cOut);
    par_ptr = libpointer('doublePtr', [radi(pore_i) rand_num(pore_i)]);    % radi dependency of reaction, R Best, variant robac for every realization, specify external parameters from dll

    % Calling BRNS
    calllib('brns','invokebrns_',cat_ptr, cbt_ptr, oc_ptr, nos_ptr, ts_ptr, fcb_ptr, rv_ptr, posx_ptr, posy_ptr, posz_ptr, p_ptr, ws_ptr, par_ptr)

    % Retrieving Data
    returnValue=get(rv_ptr, 'Value');
    if returnValue~=0, [returnValue pore_i],end
    % CHECK returnValue:
    % 0: Computation ok
    % 1: Negative concentrations occured
    % 2: Exceeding maximum newton iteration
    % 3: 1 and 2
        
    cOut (1:spec_no)    = get(oc_ptr, 'Value');
    bac(pore_i)         = round(cOut(1));        
    cNew(pore_i)        = cOut(2);
            
end

unloadlibrary('brns')                                                      % close BRNS library

spec = struct('chem', cNew, 'bac', bac);

end