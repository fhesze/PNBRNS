function TransportData = getTransport(FlowData, GeometryData, TransportCoeffs)

%%-- V 1.6.1 --------------------------------------------------------------
% Advection, Martin Explicite
% we consider only pores that feed the specific pore

%%** NOTE for change of BC ************************************************
% after changing the boundary condition and flow from Left to Right, we
% don't encounter the parallel pores across the fluid flow. so some parts
% like pore_adv removed...

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n           = GeometryData.GeometryCoeffs.xIncr;
m           = GeometryData.GeometryCoeffs.yIncr;
pore_len    = GeometryData.GeometryCoeffs.LengthOfPore;
node_no     = GeometryData.NodeData.NumberOfNodes;
nodes       = GeometryData.NodeData.NodesOfPores;
pore_info   = GeometryData.PoreData.PoreInfo;
poreXY      = GeometryData.PoreData.PorePos;

flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores
pore_CS     = pi.*radi.^2;                                                 % pore cross section
pore_vol    = pi.*radi.^2.*pore_len;                                       % volume of pores

%%-- Transport parameters -------------------------------------------------

c.spec_no   = TransportCoeffs.c_spec_no;                                   % number of chemical species 
b.spec_no   = TransportCoeffs.b_spec_no;                                   % number of bacterial species
spec_no     = c.spec_no + b.spec_no;                                       % total number of species

%%-- Time parameters ------------------------------------------------------

t_end       = TransportCoeffs.t_end;
t_0         = TransportCoeffs.t_0;
dt          = TransportCoeffs.dt;
% t_i      = 0;
t_save      = 20;                                                          % 20 % dtXt_save=1

%%-- Initial and boundary conditions of the chemical species --------------

c.T0         = zeros(pores_no, c.spec_no);
% dc.X0(1)_dt  = -0.1;                                                      % BC changing by time

c.X0(1:c.spec_no) = 1.55e-6;                                                 % c.In_rand(i_radi);

% for BC Diffusion (from 1 side to other side)
nodeX0      = find(nodes(:,7)==1);
% nodeXL      = find(nodes(:,7)==2);

%%--  BRNS file fixed params in lib-file ----------------------------------% this Av should be similar and somehow equal for hetrogenous networks of normally distributed values of radi

param       = get_param;

Av_pore     = 2*2./radi;                                                   % increase artificialy the Av of 160mic pore network by 2 to make it equal to Av from exp
Area_tot    = 2*pi*sum(radi)*pore_len;
Vol_tot     = sum(pore_vol);
Av_model    = 2*Area_tot/Vol_tot;                                          % x2 for rough microbe surface
robac       = param.massbac/param.surfexp^2/param.volexp/Av_model;         % To make a new bac density for hetrogenos medium in order to have similar conditions to homogenous case

qMax        = Av_pore.*param.mumax.*robac./1000;                           % Mole/sec
thiele      = (qMax.*radi)./(param.diff*param.Km.*Av_pore);

%%-- Advection parameters -------------------------------------------------

qAdv        = abs(flux);
qAdvEff     = get_vEff(thiele).*qAdv;

%%-- Diffusion parameters -------------------------------------------------
                          
D(1:c.spec_no)  = param.diff;      

%%-- predefined variables (for efficiency) --------------------------------

% cNew        = c.T0; 
c.Old        = c.T0;  
c.New        = zeros(pores_no, c.spec_no); 
c.Adv        = zeros(pores_no, c.spec_no);                                  % 
c.Diff       = zeros(pores_no, c.spec_no); 
qDiff       = zeros(pores_no, c.spec_no); 
c.In         = zeros(pores_no, c.spec_no);                                  % inlet pore network concentration
c.Node       = zeros(node_no, c.spec_no);     
c.Save       = zeros(pores_no, t_end/dt/t_save, c.spec_no);              % concentration profiles for saving and plotting

%%=========================================================================
%%-- Central solver -------------------------------------------------------
%%=========================================================================

t        = t_0;
t_i      = 0;

while t<t_end
    
    %%-- Transport of chemical species (advection) ------------------------
    
    c.In         = get_cIn(GeometryData, FlowData, c, qAdvEff);
    
    QAdvEff = repmat(qAdvEff,1,c.spec_no);
    Pore_Vol= repmat(pore_vol,1,c.spec_no);
    c.Adv    = dt.*QAdvEff./Pore_Vol.*(c.In - c.Old);
    
    %%-- Advection (/end) -------------------------------------------------
        
    % c.Old       = c.Old + c.Adv;
    
    %%-- Transport of chemical species (diffusion) ------------------------
	   
    c.Node       = get_cNode(GeometryData, FlowData, c);
    
    c.Node(nodeX0,:) = repmat(c.X0,size(nodeX0),1);
    DD              = repmat(D,pores_no,1);

    qDiff   = c.Node(pore_info(:,1),:) - 2.*c.Old + c.Node(pore_info(:,2),:);
    c.Diff   = 2.*DD.*dt./pore_len.^2.*qDiff;
    
	%%-- Diffusion (/end) -------------------------------------------------
    
    c.New = c.Old + c.Diff + c.Adv;
    
    %%-- Reaction (BRNS 1.8) ----------------------------------------------
     
    % c.New = getReaction(FlowData, GeometryData, TransportCoeffs, c);
    % c.New = getReaction_2(FlowData, GeometryData, TransportCoeffs, c);
%     c.New = getReaction_3(FlowData, GeometryData, TransportCoeffs, c.New);
    
    %%-- Reaction (/ends) -------------------------------------------------

    %%-- Transport of microbial species -----------------------------------

    % 
    
    %%-- Transport (/ends) ------------------------------------------------   
       
    %%-- Saving -----------------------------------------------------------
    
    t   = t + dt;
    t_i = t_i + 1;
    
    if mod(t_i,t_save)==0
        c.Save(:, t_i/t_save, 1:c.spec_no) = c.New(:, 1:c.spec_no);
    end
    
    %%-- Saving (/ends) ---------------------------------------------------
    
    c.Old = c.New;                                                         % last line, passing the concentration values to the next iteration 
    
end

%%=========================================================================
%%-- Postprocessing -------------------------------------------------------
%%=========================================================================
   








%%=========================================================================
%%-- Plotting -------------------------------------------------------------
%%=========================================================================

% for t_i = 4:-1:1
%     
%     figure;
%     scatter(poreXY(:,1), poreXY(:,2), 50, c(:,end - 10*t_i,1), '.');
%     axis equal;
%     
% end

% for spec_i = 1:c.spec_no
%     
%     figure;
%     scatter(poreXY(:,1), poreXY(:,2), 50, c.Save(:,end,spec_i), '.');
%     axis equal;
%     
% end

%%-------------------------------------------------------------------------

TransportData = struct('value', c.Save);

end

%%=========================================================================
%%-- Auxillary functions --------------------------------------------------
%%=========================================================================

function cIn = get_cIn(GeometryData, FlowData, c, qAdvEff)

%%-- Geometry parameters --------------------------------------------------

n           = GeometryData.GeometryCoeffs.xIncr;
m           = GeometryData.GeometryCoeffs.yIncr;
nodes       = GeometryData.NodeData.NodesOfPores;
pore_info   = GeometryData.PoreData.PoreInfo;
poreXY      = GeometryData.PoreData.PorePos;

%%-- Flow parameters ------------------------------------------------------

flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores

%%-- Transport parameters -------------------------------------------------

cOld        = c.Old;
cX0         = c.X0;
c_spec_no   = c.spec_no;                                                    % number of chemical species 
                                          
for pore_i = 1:pores_no

    sign_flux   = sign(flux(pore_i));
    
    if sign_flux == 1                                                      % Type I influx
        node_edge   = pore_info(pore_i,1);
    elseif sign_flux == -1                                                 % Type II influx
        node_edge    = pore_info(pore_i,2);
    end
                    
    if  pore_info(pore_i,4) == 1                                           % node at the inlet of the geometry
        
        cIn(pore_i,1:c_spec_no) = cX0; 
        
    elseif pore_info(pore_i,4) ~= 1                                        % node at the inner geometry
                   
        if sign_flux                                             
                
            pore_in = nodes(node_edge, [1,2,3]);
                
            pore_in(pore_in == 0)                       = [];
            pore_in(pore_in == pore_i)                  = [];
            pore_in(sign(flux(pore_in)) == sign_flux)   = [];
                
            for spec_i = 1:c_spec_no
                cIn(pore_i,spec_i) = sum(cOld(pore_in,spec_i).*qAdvEff(pore_in))/sum(qAdvEff(pore_in));                    
            end
                           
        else                                                               % no influx          
             cIn(pore_i,1:c_spec_no) = 0; 
        end
        
    end

end

end

function cNode = get_cNode(GeometryData, FlowData, c)

nodes       = GeometryData.NodeData.NodesOfPores;
node_no     = GeometryData.NodeData.NumberOfNodes;
pore_CS     = pi.*FlowData.radii.^2;                                       % pore cross section

cNode       = zeros(node_no, c.spec_no);

for node_i = 1:node_no
        
    neighb_pores            = nodes(node_i,[1,2,3]);
    cut_pores               = neighb_pores<=0;
    neighb_pores(cut_pores) = [];
        
    c1                      = pore_CS(neighb_pores);
    C1                      = repmat(c1, 1, c.spec_no);
    c2                      = sum(pore_CS(neighb_pores));
           
    cNode(node_i,:)         = sum(c.Old(neighb_pores,1:end).*C1)./c2;
       
end

end

function param = get_param

fid                 = fopen('params/params.txt');

foo         = fgetl(fid);
foo(1:5)    = [];
param.Km    = str2num(foo);

foo         = fgetl(fid);
foo(1:8)    = [];
param.mumax = str2num(foo);

foo          = fgetl(fid);
foo(1:7)     = [];
param.diff  = str2num(foo);

foo       = fgetl(fid);
foo(1:10) = [];
param.massbac  = str2num(foo);

foo       = fgetl(fid);
foo(1:10) = [];
param.surfexp = str2num(foo);

foo        = fgetl(fid);
foo(1:9)   = [];
param.volexp = str2num(foo);

fclose(fid);

end

function vEff = get_vEff(thiele)

vEff = 1 + (0.42./(exp((-log(thiele) + log(.02))/.95).*100 + 1));          % Effective velocity results from reaction in a single pores

end
