function FlowData = getFlow(NetworkData, GeometryData, FlowCoeffs)

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n               = GeometryData.GeometryCoeffs.xIncr;
m               = GeometryData.GeometryCoeffs.yIncr;
porelength      = GeometryData.GeometryCoeffs.LengthOfPore;
node_no         = GeometryData.NodeData.NumberOfNodes;
node_realXY     = GeometryData.NodeData.Nodes;
x_actnode       = GeometryData.NodeData.XNode;
y_actnode       = GeometryData.NodeData.YNode;
node_pores_nodes= GeometryData.NodeData.NodesOfPores;
pore_pos        = GeometryData.PoreData.PorePosition;
pore_realXY     = GeometryData.PoreData.Pores;

radi            = NetworkData.value;
desired_DH      = FlowCoeffs.delta_p;                                      % noch "Uberdenken

% load('radi_sigma2e-9_lx005','radi_het','desired_DH');
% radi            = radi_het(:,1);
% desired_DH(2:5) = [];

pores_no        = n/2*(m-1) + (n/2-1)*fix(m/2);    
xpore           = zeros(pores_no,1);
ypore           = zeros(pores_no,1);
logx            = x_actnode(pore_pos(:,2))-x_actnode(pore_pos(:,1));
xpore(logx==-1) = x_actnode(pore_pos((logx==-1),2))+.5;
xpore(logx==1)  = x_actnode(pore_pos((logx==1),1))+.5;

mio_water       = 8.94e-4;
Cond            = pi.*radi.^4./(8*mio_water*porelength);
vol             = pi.*radi.^2.*porelength;
pore_area       = pi.*radi.^2;

nodes_notBC     = node_no-m;                                               % only for even m

%%=========================================================================
%%-- Zentraler L"oser -----------------------------------------------------
%%=========================================================================

downdiag = zeros(nodes_notBC,1);
downl = downdiag;
downr = downdiag;

updiag = downdiag;
upr = updiag;
upl = updiag;

bvector  = zeros(nodes_notBC,1);
maindiag = zeros(nodes_notBC,1);

iup = m;
idown = 0;
idiag = 0;

for inod = m/2+1:node_no-m/2
    
    iup   = iup+1;
    idown = idown+1;
    idiag = idiag+1;
    
    if rem(x_actnode(inod),2) == 0 
        updiag(idiag) = -Cond(node_pores_nodes(inod,3));
        % updiag(idiag)=node_pores_nodes(inod,3+3);                        % to check if its right node
       
        if idown == node_pores_nodes(inod,4)
            downdiag(idiag) = -Cond(node_pores_nodes(inod,1));
            if node_pores_nodes(inod,2)~=0,downr(idiag) = -Cond(node_pores_nodes(inod,2));end
            % downr(idiag)=node_pores_nodes(inod,2+3);                      % to check if its right node
            % downdiag(idiag)=node_pores_nodes(inod,1+3);
        else
            downdiag(idiag)=-Cond(node_pores_nodes(inod,2));
            if node_pores_nodes(inod,1)~=0,downl(idiag) = -Cond(node_pores_nodes(inod,1));end
            % downl(idiag)=node_pores_nodes(inod,1+3);                      % to check if its right node
            % downdiag(idiag)=node_pores_nodes(inod,2+3);
        end
        
    else
        downdiag(idiag)=-Cond(node_pores_nodes(inod,2));
        % downdiag(idiag)=node_pores_nodes(inod,2+3);                       % to check if its right node
        if iup==node_pores_nodes(inod,4)
            updiag(idiag)=-Cond(node_pores_nodes(inod,1));
            if node_pores_nodes(inod,3)~=0,upr(idiag) = -Cond(node_pores_nodes(inod,3));end
            % upr(idiag)=node_pores_nodes(inod,3+3);                        % to check if its right node
            % updiag(idiag)=node_pores_nodes(inod,1+3);
        else
            updiag(idiag) = -Cond(node_pores_nodes(inod,3));
            if node_pores_nodes(inod,1)~=0,upl(idiag) = -Cond(node_pores_nodes(inod,1));end
            % upl(idiag)=node_pores_nodes(inod,1+3);%to check if its right node
            % updiag(idiag)=node_pores_nodes(inod,3+3);
        end
    end
end

% BC: the head pressure values in starting and ending nodes
% each applied to the first n/2 and last n/2 nodes in the structure
head_begin = desired_DH;%62.15 for 160mic%248.6 for 80mic;%deltaH(i_radi);
head_end   = 0;

% downl is always zero for bv_begin and so is upl for bv_end, so we removed them in next lines but keeping the original command
% bv_end   = [upl(end-m/2+1:end),updiag(end-m/2+1:end),upr(end-m/2+1:end)];
% bv_begin = [downl(1:m/2),downdiag(1:m/2),downr(1:m/2)];
bv_end   = [updiag(end-m/2+1:end),upr(end-m/2+1:end)];
bv_begin = [downdiag(1:m/2),downr(1:m/2)];
        
for inod = 1:m/2
    
    bvb_nozero = bv_begin(inod,:);
    bvn_nozero = bv_end(end-inod+1,:);
    % bvb_nozero(bvb_nozero==0) = [];
    % bvn_nozero(bvn_nozero==0) = [];
    
    % minus sign (-) because bvector is on the other side of the equation
    bvector(inod) = -sum(bvb_nozero)*head_begin;
    bvector(end-inod+1) = -sum(bvn_nozero)*head_end;
            
end

indmaindiag = node_pores_nodes(m/2+1:end-m/2,[1 2 3]);

for inod = 1:nodes_notBC
            
    none_zero_pores                     = indmaindiag(inod,:);
    none_zero_pores(none_zero_pores==0) = [];
    maindiag(inod,1)                    = sum(Cond(none_zero_pores));
            
end
%%-- Matrix solved --------------------------------------------------------

Amat3 = spdiags([upr updiag upl maindiag downr downdiag downl], [-(m/2+1) -m/2 -(m/2-1) 0 m/2-1 m/2 m/2+1], nodes_notBC, nodes_notBC);
% Amat = spdiags([downl downdiag downr maindiag upl updiag upr], [-(m/2+1) -m/2 -(m/2-1) 0 m/2-1 m/2 m/2+1], nodes_notBC, nodes_notBC);
head_calc = Amat3\bvector;

%%-- Head values over medium & Fluxes -------------------------------------

head_BC1 = repmat(head_begin,m/2,1);
head_BC2 = repmat(head_end,m/2,1);
head     = [head_BC1; head_calc; head_BC2];

%flowartes inside the pores, predefind flow directions from pores_pos
flux = -Cond.*(head(pore_pos(:,2)) - head(pore_pos(:,1)));

%%=========================================================================
%%-- Plotting -------------------------------------------------------------
%%=========================================================================

% figure;
% scatter(node_realXY(:,1), node_realXY(:,2), 50, head, '.');
% axis equal;
figure;
scatter(pore_realXY(:,1), pore_realXY(:,2), 50, flux, '.');
axis equal;

% plotting heads over the domain
% z = reshape(head, m/2, length(head)/(m/2));
% y = reshape(y_actnode, m/2, length(y_actnode)/(m/2));
% x = reshape(x_actnode, m/2, length(x_actnode)/(m/2));
 
% figure;
% surf(x,y,z)                                                                % heads
% contour(x,y,z);                                                            % heads in nodes

% figure;
% plot3k_lable({xpore ypore abs(flux)},[],[],{},6);                          % fluxes in pores

% figure;
% plot3k_lable({xpore ypore Cond},[],[],{},6);                               % Conductivity of pores

%%=========================================================================
%%-- Postprocessing -------------------------------------------------------
%%=========================================================================

%some commands to find out the total inlet and oulet Q
poreX0          = find(pore_pos(:,4)==1);
poreXL          = find(pore_pos(:,4)==2);
inletQ          = sum(abs(flux(poreX0)));
outletQ         = sum(abs(flux(poreXL)));

vol_tot         = sum(vol);
desired_Q       = .014*vol_tot;

desired_DH      = desired_Q*head_begin/inletQ;
flux_adv_abs    = abs(flux);

%%=========================================================================
%%-- Saving -------------------------------------------------------------
%%=========================================================================

%save(radi_file,'desired_DH','flux_adv_abs','-append')
%save(radi_file,'flux_adv_abs','-append')

FlowData = struct('value', flux, 'radii', radi);

end
