function cNew = getReaction(FlowData, GeometryData, TransportCoeffs, c)

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n           = GeometryData.GeometryCoeffs.xIncr;
m           = GeometryData.GeometryCoeffs.yIncr;
pore_len    = GeometryData.GeometryCoeffs.LengthOfPore;
node_no     = GeometryData.NodeData.NumberOfNodes;
nodes       = GeometryData.NodeData.NodesOfPores;
pore_pos    = GeometryData.PoreData.PoreInfo;
poreXY      = GeometryData.PoreData.PorePos;

% flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores   
pore_CS     = pi.*radi.^2;                                                 % pore cross section
pore_vol    = pi.*radi.^2.*pore_len; 

%%-- Transport Parameters -------------------------------------------------

c_spec_no   = TransportCoeffs.c_spec_no;                                   % number of chemical species 
b_spec_no   = TransportCoeffs.b_spec_no;                                   % number of bacterial species
spec_no     = c_spec_no + b_spec_no;                                       % total number of species

%%-- Time Parameters ------------------------------------------------------

t_end       = TransportCoeffs.t_end;
t_begin     = TransportCoeffs.t_0;
dt          = TransportCoeffs.dt;
% time_i      = 0;
% time_save   = 20; 

%%-------------------------------------------------------------------------

cOld        = c;
% bac         = spec.bac;

%%--  BRNS file fixed params in lib-file ----------------------------------% this Av should be similar and somehow equal for hetrogenous networks of normally distributed values of radi

Km              = 231e-9;                                                  % Mole
mumax           = 0.326e-9;                                                % mole/mgBac/sec
Diff_BRNS       = 6e-10;                                                   % equal to diffusivity coffs of species
massbac         = 0.248;                                                   % massbac_vec(i_radi);%for Sensitive analysis %mg
surfexp         = 606.6e-4;                                                % m^2
volexp          = 2.45e-6;                                                 % m^3

Area_tot        = 2*pi*sum(radi)*pore_len;
Vol_tot         = sum(pore_vol);
Av_model        = 2*Area_tot/Vol_tot;                                      % x2 for rough microbe surface
robac           = massbac/surfexp*surfexp/volexp/Av_model;                 % To make a new bac density for hetrogenos medium in order to have similar conditions to homogenous case

% Av_pore         = 2*2./radi;                                               % increase artificialy the Av of 160mic pore network by 2 to make it equal to Av from exp
% qMax            = Av_pore.*mumax.*robac./1000;                             % Mole/sec
% thiele          = (qMax.*radi)./(Diff_BRNS*Km.*Av_pore);
% Veff            = 1+(0.42./(exp((-log(thiele) + log(.02))/.95)*100 + 1));  % Effective velocity results from reaction in a single pores

% qAdvEff         = Veff.*qAdv;
% mio_water       = 8.94e-4;                                                 % water viscosity
% Cond            = pi.*radi.^4./(8*mio_water*pore_len);                     % conductivity profile based on Haggen-P equation

%%-------------------------------------------------------------------------

% pos_x           = 0;
% pos_y           = 0;
% pos_z           = 0;
% posx_ptr        = libpointer('doublePtr', pos_x);
% posy_ptr        = libpointer('doublePtr', pos_y);
% posz_ptr        = libpointer('doublePtr', pos_z);
% 
% fixedConcBound  = zeros(spec_no);                                          % must have entry for every species!
% porosity        = 0;
% waterSat        = 0;
% fcb_ptr         = libpointer('int32Ptr', fixedConcBound);
% p_ptr           = libpointer('doublePtr', porosity);
% ws_ptr          = libpointer('doublePtr', waterSat);
% 
% returnValue     = 0;                                                       % -1;
% rv_ptr          = libpointer('int32Ptr', returnValue);
% nos_ptr         = libpointer('int32Ptr', c_spec_no);
% 
% timeStep        = dt;                                                      % for now we take same timestep from advection for reaction
% ts_ptr          = libpointer('doublePtr', timeStep);


%%-------------------------------------------------------------------------

% rand_num        = rand(pores_no,1);

%%-------------------------------------------------------------------------

% fid             = fopen('arch/conf.txt');
% os              = fgetl(fid);
% os(1:5)         = [];
% file            = fgetl(fid);
% file(1:7)       = [];
% header          = fgetl(fid);
% header(1:9)     = [];
% fclose(fid);
% 
% path            = getLibPath(os, 'Scenario II');
% lib_file        = [path file];
% lib_header      = [path header];
% 
% [notfnd, warn]  = loadlibrary(lib_file, lib_header);

%%-------------------------------------------------------------------------

for pore_i=1:pores_no
        
%         cAfterTransport(1:c_spec_no)= cNew(pore_i,1:c_spec_no);
%         cInit                       = cAfterTransport;                     % used as initial guess for newton iteration
%         cOut(1:c_spec_no)           = 0.0; 
%         rrVec(1:c_spec_no)          = 0.0; 
%        
%         cat_ptr = libpointer('doublePtr', cAfterTransport);
%         cit_ptr = libpointer('doublePtr', cInit);
%         co_ptr  = libpointer('doublePtr', cOut);
%         rrv_ptr = libpointer('doublePtr', rrVec);
%         par_ptr = libpointer('doublePtr', [radi(pore_i) robac 1.5*Km]);    % radi dependency of reaction, R Best, variant robac for every realization
%         
%         % Calling BRNS
%         calllib('brns','invokebrns_',cat_ptr, cit_ptr, co_ptr, nos_ptr, ts_ptr, fcb_ptr, rv_ptr, posx_ptr, posy_ptr, posz_ptr, p_ptr, ws_ptr, par_ptr);
% 
%         % Retrieving Data
%         returnValue = get(rv_ptr, 'Value');
%         if returnValue~=0, [returnValue pore_i],end
%         % CHECK returnValue:
%         % 0: Computation ok
%         % 1: Negative concentrations occured
%         % 2: Exceeding maximum newton iteration
%         % 3: 1 and 2
% 
%         cOut            = get(co_ptr, 'Value');
%         rrVec           = get(rrv_ptr, 'Value');
        
        cDelta = cOld(pore_i,:)*mumax;
        cNew(pore_i,:)  = cOld(pore_i,:) - cDelta;
        
end

% unloadlibrary('brns')                                                      % close BRNS library

% spec = struct('chem', cNew);

end
