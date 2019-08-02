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
nodeXY      = GeometryData.NodeData.Nodes;
pore_info   = GeometryData.PoreData.PoreInfo;
poreXY      = GeometryData.PoreData.PorePos;

flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores
pore_CS     = pi.*radi.^2;                                                 % pore cross section
pore_vol    = pi.*radi.^2.*pore_len;                                       % volume of pores

c.spec_no   = TransportCoeffs.c_spec_no;                                   % number of chemical species 
b.spec_no   = TransportCoeffs.b_spec_no;                                   % number of bacterial species
spec_no     = c.spec_no + b.spec_no;                                       % total number of species

%%-- Time parameters ------------------------------------------------------

t_end       = TransportCoeffs.t_end;
t_0         = TransportCoeffs.t_0;
dt          = TransportCoeffs.dt;
t           = t_0:dt:t_end;
t_save      = 20;                                                          

%%-- Initial and boundary conditions of the chemical species --------------

c.T0         = zeros(pores_no, c.spec_no);                                  % inital condition
c.X0(1:c.spec_no) = 1.55e-6;                                                % boundary condition

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

c.Old       = c.T0;                                                        % concentration on the old time step
c.OldDown   = c.T0;                         

c.New       = zeros(pores_no, c.spec_no);                                  
c.NewDown   = zeros(pores_no, c.spec_no);

c.Adv       = zeros(pores_no, c.spec_no);                                  
c.AdvDown   = zeros(pores_no, c.spec_no);

c.Diff      = zeros(pores_no, c.spec_no); 
c.DiffDown  = zeros(pores_no, c.spec_no);

qDiff       = zeros(pores_no, c.spec_no); 
qDiffDown   = zeros(pores_no, c.spec_no);

c.In        = zeros(pores_no, c.spec_no);                                  % inlet pore network concentration
c.InDown    = zeros(pores_no, c.spec_no);

c.Node      = zeros(node_no, c.spec_no);     
c.NodeDown  = zeros(node_no, c.spec_no);

c.Save      = zeros(pores_no, t_end/dt/t_save, c.spec_no);                 % concentration profiles for saving and plotting

%%=========================================================================
%%-- Central solver -------------------------------------------------------
%%=========================================================================
   
PoreIn = get_poreIn(GeometryData, FlowData);

for t_i = t
    
    %%-- Transport of chemical species (advection) ------------------------
    
    QAdvEff     = repmat(qAdvEff, 1, c.spec_no);
    Pore_Vol    = repmat(pore_vol, 1, c.spec_no);
    
    c.In        = get_cIn(c, qAdvEff, PoreIn, 'upper');
    c.InDown    = get_cIn(c, qAdvEff, PoreIn, 'lower');
    
    c.Adv       = dt.*QAdvEff./Pore_Vol.*(c.In - c.Old);
    c.AdvDown   = dt.*QAdvEff./Pore_Vol.*(c.InDown - c.OldDown);
    
    %%-- Advection (/end) -------------------------------------------------
    
    % c.Old       = c.Old + c.Adv;
%     c.OldDown  = c.OldDown + c.AdvDown;
    
    %%-- Transport of chemical species (longitudinal diffusion) -----------
	   
    DD          = repmat(D,pores_no,1);

    c.Node       = get_cNode(GeometryData, FlowData, c, 'upper');
    c.NodeDown   = get_cNode(GeometryData, FlowData, c, 'lower');
    
    c.Node(nodeX0,:) = repmat(c.X0,size(nodeX0),1);
    c.NodeDown(nodeX0,:) = repmat(c.X0,size(nodeX0),1);

    qDiff       = c.Node(pore_info(:,1),:) - 2.*c.Old + c.Node(pore_info(:,2),:);
    c.Diff       = 2.*DD.*dt./pore_len.^2.*qDiff;
    
    qDiffDown   = c.NodeDown(pore_info(:,1),:) - 2.*c.OldDown + c.NodeDown(pore_info(:,2),:);
    c.DiffDown   = 2.*DD.*dt./pore_len.^2.*qDiffDown;
    
	%%-- Longitudinal Diffusion (/end) ------------------------------------
    
    %%-- Transport of chemical species ( transversal diffusion) -----------
    
%     k_tr = pi^2/4.*D.*Av_pore./radi;
% 
%     for pore_i = 1:pores_no
%         
%         dy      = radi(pore_i);
%         qTrDiff = pi.^2./2.*D.*(c.Old(pore_i,:) - c.OldDown(pore_i,:))./dy;
% 
%         c.TrDiff(pore_i,:)    = - dt./dy.*qTrDiff;
%         c.TrDiffDown(pore_i,:)   =   dt./dy.*qTrDiff;
%         
%     end
    
    %%-- Transversal Diffusion (/end) -------------------------------------
     
    c.New       = c.Old + c.Diff + c.Adv; % + c.TrDiff;
    c.NewDown   = c.OldDown + c.DiffDown + c.AdvDown; % + c.TrDiffDown;
    
    %%-- Reaction (BRNS 1.8) ----------------------------------------------
           
%     c.New = getReaction(FlowData, GeometryData, TransportCoeffs, c);
%     c.New = getReaction_2(FlowData, GeometryData, TransportCoeffs, c);
    
%     c.New = getReaction_3(FlowData, GeometryData, TransportCoeffs, c.New);
%     c.NewDown = getReaction_3(FlowData, GeometryData, TransportCoeffs, c.NewDown);
    
    %%-- Reaction (/ends) -------------------------------------------------

    %%-- Transport of microbial species -----------------------------------

    % 
    
    %%-- Transport (/ends) ------------------------------------------------   
       
    %%-- Saving -----------------------------------------------------------
    
    t_index = t_i/dt + 1;

    if mod(t_index, t_save)==0 
        c.Save(:, t_index/t_save, 1:c.spec_no) = c.NewDown(:, 1:c.spec_no);
    end
    
    %%-- Saving (/ends) ---------------------------------------------------
    
    c.Old = c.New;                                                         % last line, passing the concentration values to the next iteration 
    c.OldDown = c.NewDown; 
    
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

% figure;
% scatter(poreXY(:,1), poreXY(:,2), 200, FlowData.value, 'filled');
% axis equal;

figure;
scatter(poreXY(:,1), poreXY(:,2), 200, PoreIn.flag, 'filled');
axis equal;

% figure;
% scatter(poreXY(:,1), poreXY(:,2), 200,  pore_info(:, 3), 'filled');
% axis equal;

% figure(2);
% hold on;
% for pore_i = 1:pores_no
%     
%     x_pos = poreXY(pore_i,1);
%     y_pos = poreXY(pore_i,2);
%     
%     if  strcmp(PoreIn.type(pore_i), 'inlet_pore')
%         color = 'b';
%     elseif strcmp(PoreIn.type(pore_i), 'boundary_pore')
%         color = 'm';
%     elseif strcmp(PoreIn.type(pore_i), 'tilted_pore_double_inlet')
%         color = 'g';
%     elseif strcmp(PoreIn.type(pore_i), 'straight_pore_double_inlet')
%         color = 'y';
%     elseif strcmp(PoreIn.type(pore_i), 'tilted_pore_single_inlet')
%         color = 'r';
%     end
%         
%     if (pore_info(pore_i, 3) == 1)
%         plot(x_pos, y_pos, 'v', 'MarkerFaceColor', color, 'MarkerSize',20);
%     elseif  (pore_info(pore_i, 3) == 2)
%         plot(x_pos, y_pos, '^', 'MarkerFaceColor', color, 'MarkerSize',20);
%     elseif (pore_info(pore_i, 3) == 3)
%         plot(x_pos, y_pos, 's', 'MarkerFaceColor', color, 'MarkerSize',20);
%     end
% 
% end
% hold off;
% axis equal;

%%-------------------------------------------------------------------------

TransportData = struct('value', c.Save);

end

%%=========================================================================
%%-- Auxillary functions --------------------------------------------------
%%=========================================================================

function PoreIn = get_poreIn(GeometryData, FlowData)

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
qAdv        = abs(flux);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores

%%-- Transport parameters -------------------------------------------------

for pore_i = 1:pores_no

    sign_flux   = sign(flux(pore_i));
    PoreIn.pore_info(pore_i,:) = pore_info(pore_i,:);
    
    % Type I influx
    if sign_flux == 1                                                      
        node_edge = pore_info(pore_i, 1);
    % Type II influx    
    elseif sign_flux == -1                                                 
        node_edge = pore_info(pore_i, 2);
    end
       
    % node at the inlet of the geometry
    if  pore_info(pore_i,4) == 1
        PoreIn.type{pore_i,:} = 'inlet_pore';
        PoreIn.node_type{pore_i,:} = 'type_1';
        PoreIn.flag(pore_i) = -1;
        PoreIn.data{pore_i,:} = -1;
        PoreIn.pore_pos_y(pore_i) = poreXY(pore_i,2);
        
    % node at the inner geometry
    elseif pore_info(pore_i,4) ~= 1                                        
        
        % flux is bigger than zero
        if sign_flux                                             
                
            pore_in = nodes(node_edge, [1,2,3]);
            
            pore_in(pore_in == 0)                       = [];
            pore_in(pore_in == pore_i)                  = [];
            pore_in(sign(flux(pore_in)) == sign_flux)   = [];
            
%             if (sign_flux == 1) && (pore_info(pore_i,3) == 1)
%                 PoreIn.type{pore_i,:} = 'tilted_pore_double_inlet';
%                 PoreIn.flag(pore_i) = 0.9;
%             elseif (sign_flux == 1) && (pore_info(pore_i,3) == 2)
%                 PoreIn.type{pore_i,:} = 'tilted_pore_double_inlet';
%                 PoreIn.flag(pore_i) = 1.1;  
%             end
            
            
            if length(find(nodes(node_edge, [1,2,3]))) == 2
                PoreIn.type{pore_i,:} = 'boundary_pore';
                PoreIn.node_type{pore_i,:} = 'type_2';
                PoreIn.flag(pore_i) = 4.1;
                
            elseif (length(pore_in) == 1) && (pore_info(pore_i,3) == 1)              
                PoreIn.type{pore_i,:} = 'uptilted_pore_single_inlet';
                PoreIn.flag(pore_i) = 2.9;
                
                % type III node
                if pore_info(pore_in,3) == 3
                    PoreIn.node_type{pore_i,:} = 'type_3';
                    pore_info(nodes(node_edge, [1]),3);
                    q2 = qAdv( nodes(node_edge, [1]) );
                    pore_info(nodes(node_edge, [2]),3);
                    q3 = qAdv( nodes(node_edge, [2]) );
                    pore_info(nodes(node_edge, [3]),3);
                    q1 = qAdv( nodes(node_edge, [3]) );
                    
                    p = q2/(q1 + q2);
                    q = q1/(q1 + q2);
                % type VI node
                elseif pore_info(pore_in,3) == 2
                    PoreIn.node_type{pore_i,:} = 'type_6';
                    pore_info(nodes(node_edge, [1]),3);
                    q1 = qAdv( nodes(node_edge, [1]) );
                    pore_info(nodes(node_edge, [2]),3);
                    q2 = qAdv( nodes(node_edge, [2]) );
                    pore_info(nodes(node_edge, [3]),3);
                    q3 = qAdv( nodes(node_edge, [3]) );
                    
                    p = q1/(q1 + q3);
                    q = q3/(q1 + q3);
                end
                
                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
                
            elseif (length(pore_in) == 1) && (pore_info(pore_i,3) == 2)
                PoreIn.type{pore_i,:} = 'downtilted_pore_single_inlet';
                PoreIn.flag(pore_i) = 3.1;
                
                % type III node
                if pore_info(pore_in,3) == 3
                    PoreIn.node_type{pore_i,:} = 'type_3';
                    pore_info(nodes(node_edge, [1]),3);
                    q2 = qAdv( nodes(node_edge, [1]) );
                    pore_info(nodes(node_edge, [2]),3);
                    q3 = qAdv( nodes(node_edge, [2]) );
                    pore_info(nodes(node_edge, [3]),3);
                    q1 = qAdv( nodes(node_edge, [3]) );
                
                    p = q2/(q1 + q2);
                    q = q1/(q1 + q2);
                % type VIII node
                elseif pore_info(pore_in,3) == 1
                    PoreIn.node_type{pore_i,:} = 'type_8';
                    pore_info(nodes(node_edge, [1]),3);
                    q1 = qAdv( nodes(node_edge, [1]) );
                    pore_info(nodes(node_edge, [2]),3);
                    q2 = qAdv( nodes(node_edge, [2]) );
                    pore_info(nodes(node_edge, [3]),3);
                    q3 = qAdv( nodes(node_edge, [3]) );
                    
                    p = q3/(q2 + q3);
                    q = q2/(q2 + q3);
                end
                
                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
                
            elseif (sign_flux == 1) && (pore_info(pore_i,3) == 1)
                PoreIn.type{pore_i,:} = 'uptilted_pore_double_inlet';
                PoreIn.node_type{pore_i,:} = 'type_7';
                PoreIn.flag(pore_i) = 0.9;
                PoreIn.flux(pore_i,:) = flux(nodes(node_edge, [1,2,3]));

                q1 = qAdv(pore_in(1));
                q2 = qAdv(pore_in(2));
                                
                p = q1/(q1 + q2);
                q = q2/(q1 + q2);
                PoreIn.ratio_total(pore_i,:)= [p/2 p/2 q/2 q/2];

                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
                
            elseif (sign_flux == 1) && (pore_info(pore_i,3) == 2)
                PoreIn.type{pore_i,:} = 'downtilted_pore_double_inlet';
                PoreIn.node_type{pore_i,:} = 'type_5';
                PoreIn.flag(pore_i) = 1.1; 
                PoreIn.flux(pore_i,:) = flux(nodes(node_edge, [1,2,3]));
                q1 = qAdv(pore_in(1));
                q2 = qAdv(pore_in(2));
                                
                p = q1/(q1 + q2);
                q = q2/(q1 + q2);
                PoreIn.ratio_total(pore_i,:)= [p/2 p/2 q/2 q/2];

                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
                
            elseif (length(pore_in) == 1) && (pore_info(pore_i,3) == 3)
                PoreIn.type{pore_i,:} = 'straight_pore_single_inlet';
                PoreIn.flag(pore_i) = 4.0;  
                
               % type VI node
               if pore_info(pore_in,3) == 2
                   PoreIn.node_type{pore_i,:} = 'type_6';
                   pore_info(nodes(node_edge, [1]),3);
                   q1 = qAdv( nodes(node_edge, [1]) );
                   pore_info(nodes(node_edge, [2]),3);
                   q2 = qAdv( nodes(node_edge, [2]) );
                   pore_info(nodes(node_edge, [3]),3);
                   q3 = qAdv( nodes(node_edge, [3]) );
                    
                   p = q1/(q1 + q3);
                   q = q3/(q1 + q3);
                % type VIII node
                elseif pore_info(pore_in,3) == 1
                    PoreIn.node_type{pore_i,:} = 'type_8';
                    pore_info(nodes(node_edge, [1]),3);
                    q1 = qAdv( nodes(node_edge, [1]) );
                    pore_info(nodes(node_edge, [2]),3);
                    q2 = qAdv( nodes(node_edge, [2]) );
                    pore_info(nodes(node_edge, [3]),3);
                    q3 = qAdv( nodes(node_edge, [3]) );
                    
                    p = q3/(q2 + q3);
                    q = q2/(q2 + q3);
                end
                
                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
                
            elseif (sign_flux == -1)
                PoreIn.type{pore_i,:} = 'straight_pore_double_inlet';
                PoreIn.node_type{pore_i,:} = 'type_4';
                PoreIn.flag(pore_i) = 2.0;
                PoreIn.flux(pore_i,:) = flux(nodes(node_edge, [1,2,3]));
                q1 = qAdv(pore_in(1));
                q2 = qAdv(pore_in(2));

                p = q1/(q1 + q2);
                q = q2/(q1 + q2);
                PoreIn.ratio_total(pore_i,:)= [p/2 p/2 q/2 q/2];

                if q > p
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 (1/2 - p) 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 p/2 q/2];
                elseif p > q
                    PoreIn.ratio_down(pore_i,:) = [p/2 q/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 (1/2 - q) q/2 q/2];
                else
                    PoreIn.ratio_down(pore_i,:) = [p/2 p/2 0.0 0.0];
                    PoreIn.ratio_up(pore_i,:)   = [0.0 0.0 q/2 q/2];
                end
%                 tmp = sum(PoreIn.ratio_down(pore_i,:) +  PoreIn.ratio_up(pore_i,:)) 
            end
                        
            PoreIn.data{pore_i,:} = pore_in;
            PoreIn.p(pore_i) = p;
            PoreIn.q(pore_i) = q;
            
        % flux equals zero, i.e. no influx
        else     
            PoreIn.type{pore_i,:} = 'no_flux_pore';
            PoreIn.flag(pore_i) = 0; 
            PoreIn.data{pore_i,:} = 0;
        end        
    end
end

end

function cIn = get_cIn(c, qAdvEff, PoreIn, flag)

pore_no = length(PoreIn.flag);

for pore_i = 1:pore_no
    
    % pore at the inlet of the geometry
    if  strcmp(PoreIn.type(pore_i), 'inlet_pore')
        cIn(pore_i,:) = get_cIn_inlet_inlet(c,qAdvEff,PoreIn,pore_i,flag);
    % pore at the upper or lower boundary of the geometry
    elseif strcmp(PoreIn.type(pore_i), 'boundary_pore')
        cIn(pore_i,:) = get_cIn_boundary_inlet(c,qAdvEff,PoreIn,pore_i,flag);
    
    % internal pore with influx
    elseif strcmp(PoreIn.type(pore_i), 'uptilted_pore_double_inlet')                                           
        cIn(pore_i,:) = get_cIn_double_inlet(c,qAdvEff,PoreIn,pore_i,flag); 
    elseif strcmp(PoreIn.type(pore_i), 'downtilted_pore_double_inlet')
        cIn(pore_i,:) = get_cIn_double_inlet(c,qAdvEff,PoreIn,pore_i,flag); 
    elseif strcmp(PoreIn.type(pore_i), 'straight_pore_double_inlet')
        cIn(pore_i,:) = get_cIn_double_inlet(c,qAdvEff,PoreIn,pore_i,flag); 
    elseif strcmp(PoreIn.type(pore_i), 'uptilted_pore_single_inlet')
        cIn(pore_i,:) = get_cIn_single_inlet(c,qAdvEff,PoreIn,pore_i,flag);
    elseif strcmp(PoreIn.type(pore_i), 'downtilted_pore_single_inlet')
        cIn(pore_i,:) = get_cIn_single_inlet(c,qAdvEff,PoreIn,pore_i,flag);
    elseif strcmp(PoreIn.type(pore_i), 'straight_pore_single_inlet')
        cIn(pore_i,:) = get_cIn_single_inlet(c,qAdvEff,PoreIn,pore_i,flag);
    
    % internal pore without influx
    elseif strcmp(PoreIn.type(pore_i), 'no_flux_pore')
        cIn(pore_i,1:c.spec_no) = 0;
    else
        'no pore type specified'
    end
    
end

end

function cNode = get_cNode(GeometryData, FlowData, c, flag)

nodes       = GeometryData.NodeData.NodesOfPores;
node_no     = GeometryData.NodeData.NumberOfNodes;
pore_CS     = pi.*FlowData.radii.^2;                                       % pore cross section

if flag == 'upper'
    cOld = c.Old;
elseif flag == 'lower'
    cOld = c.OldDown;
end

cNode       = zeros(node_no, c.spec_no);

for node_i = 1:node_no
        
    neighb_pores            = nodes(node_i,[1,2,3]);
    cut_pores               = neighb_pores<=0;
    neighb_pores(cut_pores) = [];
        
    c1                      = pore_CS(neighb_pores);
    C1                      = repmat(c1, 1, c.spec_no);
    c2                      = sum(pore_CS(neighb_pores));
           
    cNode(node_i,:)         = sum(cOld(neighb_pores,1:end).*C1)./c2;
       
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

function cIn = get_cIn_inlet_inlet(c, qAdvEff, PoreIn, pore_i, flag)

pore_pos_y = PoreIn.pore_pos_y(pore_i);
cIn = zeros(1, c.spec_no);

% if pore_pos_y < 0.01 
%     cIn(1,1:c.spec_no) = 0.0;
% elseif pore_pos_y > 0.02
%     cIn(1,1:c.spec_no) = 0.0;
% else
%     cIn(1,1:c.spec_no) = c.X0;
% end

cIn(1,1:c.spec_no) = c.X0;

end

function cIn = get_cIn_boundary_inlet(c, qAdvEff, PoreIn, pore_i, flag)

pore_in = cell2mat(PoreIn.data(pore_i,:));
q3 = qAdvEff(pore_in);  
cIn = zeros(1, c.spec_no);

if flag == 'upper'
    cIn(1,:) = c.Old(pore_in(1),:).*q3/q3;
elseif flag == 'lower'
    cIn(1,:) = c.OldDown(pore_in(1),:).*q3/q3;
end 

end

function cIn = get_cIn_single_inlet(c, qAdvEff, PoreIn, pore_i, flag)

pore_in = cell2mat(PoreIn.data(pore_i,:));
c_spec_no = length(c.Old(1,:));
c_in_vec(1,:) = c.OldDown(pore_in(1),:);
c_in_vec(2,:) = c.Old(pore_in(1),:);
c_in =  (c_in_vec(1,:) + c_in_vec(2,:))./2;
% c_in(1:c_spec_no) = mean(c_in_vec)

cIn = zeros(1, c_spec_no);
q3 = qAdvEff(pore_in);  

for spec_i = 1:c_spec_no
    cIn(1,spec_i) = c_in(spec_i).*q3/q3;
end

p = PoreIn.p(pore_i);
q = PoreIn.q(pore_i);
% PoreIn.ratio_down(pore_i,:)
% p = PoreIn.ratio_down(pore_i,1)*2;
% q = PoreIn.ratio_up(pore_i,4)*2;

if strcmp(PoreIn.node_type(pore_i), 'type_3')
    if strcmp(PoreIn.type(pore_i), 'downtilted_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                c_in_1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,1);
                c_in_2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,1);
            elseif flag == 'upper'
                c_in_1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,2);
                c_in_2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,2);
            end 
            cIn(1,spec_i) = (c_in_1 + c_in_2)*2/p;
%             cIn(1,spec_i) = c_in_vec(1,spec_i);       
        end        
    elseif strcmp(PoreIn.type(pore_i), 'uptilted_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                c_in_1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,3);
                c_in_2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,3);
            elseif flag == 'upper'
                c_in_1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,4);
                c_in_2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,4);
            end
            cIn(1,spec_i) = (c_in_1 + c_in_2)*2/q;
%             cIn(1,spec_i) = c_in_vec(2,spec_i);        
        end        
    end  
elseif strcmp(PoreIn.node_type(pore_i), 'type_6')
    if strcmp(PoreIn.type(pore_i), 'uptilted_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                cIn(1,spec_i) = c_in_vec(1,spec_i);
            elseif flag == 'upper'
                tmp1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,2);
                tmp2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,2);
                cIn(1,spec_i) = (tmp1 + tmp2)*2/p;
            end        
        end        
    elseif strcmp(PoreIn.type(pore_i), 'straight_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                tmp1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,3);
                tmp2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,3);
                cIn(1,spec_i) = (tmp1 + tmp2)*2/q;
            elseif flag == 'upper'
                cIn(1,spec_i) = c_in_vec(2,spec_i);
            end        
        end        
    end 
elseif strcmp(PoreIn.node_type(pore_i), 'type_8')
    if strcmp(PoreIn.type(pore_i), 'straight_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                cIn(1,spec_i) = c_in_vec(1,spec_i);
            elseif flag == 'upper'
                tmp1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,2);
                tmp2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,2);
                cIn(1,spec_i) = (tmp1 + tmp2)*2/p;
            end        
        end        
    elseif strcmp(PoreIn.type(pore_i), 'downtilted_pore_single_inlet')
        for spec_i = 1:c_spec_no
            if flag == 'lower'
                tmp1 = c_in_vec(1,spec_i)*PoreIn.ratio_down(pore_i,3);
                tmp2 = c_in_vec(2,spec_i)*PoreIn.ratio_up(pore_i,3);
                cIn(1,spec_i) = (tmp1 + tmp2)*2/q;
            elseif flag == 'upper'
                cIn(1,spec_i) = c_in_vec(2,spec_i);
            end        
        end        
    end 
end

end

function cIn = get_cIn_double_inlet(c, qAdvEff, PoreIn, pore_i, flag)

% if flag == 'upper'
%     cOld = c.Old;
% elseif flag == 'lower'
%     cOld = c.OldDown;
% end

pore_in = cell2mat(PoreIn.data(pore_i,:));
c_spec_no = length(c.Old(1,:));
cIn = zeros(1, c_spec_no);

q1 = qAdvEff(pore_in(1));
q2 = qAdvEff(pore_in(2));
% q3 = qAdvEff(pore_i);  

% if strcmp(PoreIn.type(pore_i), 'straight_pore_double_inlet')
if 1
           
    c_in_vec(1,:) = c.OldDown(pore_in(1),:);
    c_in_vec(2,:) = c.Old(pore_in(1),:);
    c_in_vec(3,:) = c.OldDown(pore_in(2),:);
    c_in_vec(4,:) = c.Old(pore_in(2),:);
            
    c1 = (c_in_vec(1,:) + c_in_vec(2,:))/2;
    c2 = (c_in_vec(3,:) + c_in_vec(4,:))/2;
    
    for spec_i = 1:c_spec_no
     
        m_total = dot(PoreIn.ratio_total(pore_i,:), c_in_vec(:,spec_i));
        m_down = dot(PoreIn.ratio_down(pore_i,:), c_in_vec(:,spec_i));
        m_up = dot(PoreIn.ratio_up(pore_i,:), c_in_vec(:,spec_i));
        
        r = 2*m_down/m_total;
        s = 2*m_up/m_total;
%         r = 1;
%         s = 1;
        
        if  m_total == 0
            cIn(1,spec_i) = 0;
        else
            if flag == 'lower'
                cIn(1,spec_i) = r*(c1(spec_i)*q1 + c2(spec_i)*q2)/(q1+q2);
            elseif flag == 'upper'
                cIn(1,spec_i) = s*(c1(spec_i)*q1 + c2(spec_i)*q2)/(q1+q2);
            end
        end
    end
 
% elseif strcmp(PoreIn.type(pore_i), 'downtilted_pore_double_inlet')
%        
%     c_in_vec(1,:) = c.OldDown(pore_in(1),:);
%     c_in_vec(2,:) = c.Old(pore_in(1),:);
%     c_in_vec(3,:) = c.OldDown(pore_in(2),:);
%     c_in_vec(4,:) = c.Old(pore_in(2),:);
%     
%     c1 = (c_in_vec(1,:) + c_in_vec(2,:))/2;
%     c2 = (c_in_vec(3,:) + c_in_vec(4,:))/2;
%     
%     for spec_i = 1:c_spec_no
%         
%         cIn(1,spec_i) = (c1(spec_i)*q1 + c2(spec_i)*q2)/(q1 + q2);
% %         cIn(1,spec_i) = sum(cOld(pore_in,spec_i).*qAdvEff(pore_in))/sum(qAdvEff(pore_in));                    
%     end
%     
% elseif strcmp(PoreIn.type(pore_i), 'uptilted_pore_double_inlet')
%          
%     c_in_vec(1,:) = c.OldDown(pore_in(1),:);
%     c_in_vec(2,:) = c.Old(pore_in(1),:);
%     c_in_vec(3,:) = c.OldDown(pore_in(2),:);
%     c_in_vec(4,:) = c.Old(pore_in(2),:);
%     
%     c1 = (c_in_vec(1,:) + c_in_vec(2,:))/2;
%     c2 = (c_in_vec(3,:) + c_in_vec(4,:))/2;
%     
%     for spec_i = 1:c_spec_no
%         
%         cIn(1,spec_i) = (c1(spec_i)*q1 + c2(spec_i)*q2)/(q1 + q2);
% %         cIn(1,spec_i) = sum(cOld(pore_in,spec_i).*qAdvEff(pore_in))/sum(qAdvEff(pore_in));                    
%     end     
end

end

function vEff = get_vEff(thiele)

vEff = 1 + (0.42./(exp((-log(thiele) + log(.02))/.95).*100 + 1));          % Effective velocity results from reaction in a single pores

end
