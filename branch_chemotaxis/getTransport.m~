function TransportData = getTransport(FlowData, GeometryData, TransportCoeffs)

%%15 July 2010, Mehdi Gharasoo
%Final version of 2D Hexagonal pore network
% 7 march 2011, IBM complete

%%----------------------V 1.6.1 -----------------------------------
% Advection, Martin Explicite
% we consider only pores that feed the specific pore

%%*********** NOTE for change of BC ***************
%after changing the boundary condition and flow from Left to Right, we
%don't encounter the parallel pores across the fluid flow.

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
pore_vol    = pi.*radi.^2.*pore_len;                                       % volume of pores
pore_vol_tot= sum(pore_vol);

%%-- Transport Parameters -------------------------------------------------

c_spec_no   = TransportCoeffs.c_spec_no;                                   % number of chemical species 
b_spec_no   = TransportCoeffs.b_spec_no;                                   % number of bacterial species
spec_no     = c_spec_no + b_spec_no;                                       % total number of species

%%-- Time Parameters ------------------------------------------------------

t_end       = TransportCoeffs.t_end;
t_begin     = TransportCoeffs.t_0;
dt          = TransportCoeffs.dt;
% time_i      = 0;
time_save   = 20; 

%%-- Initial and boundary conditions of the chemical species --------------

cInput      = 0.394; % 1.553E-02;                                          % input conc of solute Mole
cInit       = zeros(pores_no, c_spec_no) + cInput;
% cInit         = 0 + (cInput-0).*rand(pores_no,1);
% cInit(5817)   = cInput;%center
% cInit(8503)   = cInput;%right distance 7
% cInit(4743)   = cInput;%right distance 4


% Boundry conditions
% dc.chem_bc1_dt= -.1;                                                         % BC changing by time

% nod_bc1     = find(nodes(:,7)==1);

% c.chem_bc1    = 0;%cInput;%1st spee
% c.chem_bc2    = 0;

%%-- Advection params -----------------------------------------------------

% qAdv            = abs(flux);

%%-- Diffusion params -----------------------------------------------------

% diffFlux        = zeros(pores_no,2,spec_no);                           
D(1:c_spec_no)  = 9e-10;
% diff1_coff      = 9e-10;    

%%-- Initial and boundary conditions of the bacterial species -------------

bacMax          = 10000;                                                   % a maximum of 10000 cells can be in a pore with avg radius 
cBacMax         = bacMax/(mean_radi^2*pi*pore_len);                        % maximum bacterial density in a pore with avg radius
bacInit         = fix(cBacMax*0.25.*pore_vol);                             % one fourth of maximum as inital bacterial distribution
bac             = bacInit;
% cBac            = bac./pore_vol;

%%-- Predefined variables (for efficiency) --------------------------------

cOld            = cInit;
cNode           = zeros(node_no,1);
Diff1_flux      = zeros(pores_no,2);
cSave           = zeros(pores_no, t_end/dt/time_save);                     % concentration profiles for saving and plotting
bacSave         = zeros(pores_no, t_end/dt/time_save);                     % bacteria profiles for saving and plotting 

%%-------------------------------------------------------------------------

spec = struct('chem', cInit, 'bac', bac);

%%=========================================================================
%%-- Central solver -------------------------------------------------------
%%=========================================================================

time        = t_begin;
time_i      = 0;

while time<t_end
    
%     time  
%     time_i
    
    %%-- Transport of chemical species (advection) ------------------------

    % no advection in this case
    
    %%-- Advection (/end) -------------------------------------------------

    %%-- Transport of chemical species (diffusion) ------------------------

    for node_i = 1:node_no

        neighb_pores            = nodes(node_i,[1,2,3]);
        cut_pores               = neighb_pores<=0;
        neighb_pores(cut_pores) = [];

        cNode(node_i) = sum(cOld(neighb_pores).*pore_CS(neighb_pores))/sum(pore_CS(neighb_pores));

    end

    Diff1_flux(:,1) = -D(1).*(cOld - cNode(pore_pos(:,1)))/(pore_len/2);
    Diff1_flux(:,2) = -D(1).*(cOld - cNode(pore_pos(:,2)))/(pore_len/2);

    cDiff = dt*(Diff1_flux(:,1) + Diff1_flux(:,2))/pore_len;

    %%-- Diffusion (/end) -------------------------------------------------

    spec.chem = cOld + cDiff; % + cAdv;

    %%-- Reaction (BRNS 1.8) ----------------------------------------------

    spec = getReaction(FlowData, GeometryData, TransportCoeffs, spec);
    
    %%-- Reaction (/end) --------------------------------------------------

    %%-- Transport of microbial species (chemotaxis) ----------------------

    spec = getChemotaxis(FlowData, GeometryData, TransportCoeffs, spec);
    
    %%-- Transport (/end) -------------------------------------------------   
    
    %%-- Saving -----------------------------------------------------------
    
    time        = time + dt;
    time_i      = time_i + 1;
        
    if time~=0 && mod(time_i,time_save)==0
        
        tSave               = time_i/time_save;
        cSave(:, tSave)     = spec.chem;        
        bacSave(:, tSave)   = spec.bac;   
        
    end

    %%-- Saving (/ends) ---------------------------------------------------
    
    cOld = spec.chem;                                                      % last line, passing the concentration values to the next iteration

end

%%=========================================================================
%%-- Postprocessing -------------------------------------------------------
%%=========================================================================

% y_trancent=pore_realY(7249)
% trancent=find(pore_realY==y_trancent);
% for i=1:15;conc_trancent(:,i)=cSave(trancent,i*10);end
% for i=1:15;bac_trancent(:,i)=bacSave(trancent,i*10);end
% xlswrite('chemo_4_6.xls',[pore_realX(trancent), bac_trancent, conc_trancent],'chemoattract')

% size(cSave)

%%=========================================================================
%%-- Plotting -------------------------------------------------------------
%%=========================================================================

figure;
scatter(poreXY(:,1), poreXY(:,2), 50, cSave(:,end), '.');
axis equal;
% axis([0 max(poreXY(:,1)) 0 max(poreXY(:,2))]);

%%-------------------------------------------------------------------------

TransportData = struct('value', cSave);

end
