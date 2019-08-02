function GeometryData = getGeometry(GeometryCoeffs)

%% Creating a two-dimensional hexagonal pore network
%% The size is determined by the parameters 'm', 'n' and 'pore_lenght'
%% The return values are the coordinates of the pores
 
%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n           = GeometryCoeffs.xIncr;
m           = GeometryCoeffs.yIncr;
porelength  = GeometryCoeffs.LengthOfPore;

nodes       = zeros(m,n);                                                  % in the matrix, m represents y-axis and n does x-axis
actnodes    = m*n/2;                                                       % only when n is factor of 4

x_actnode   = zeros(actnodes,1);
y_actnode   = zeros(actnodes,1);

yx_actnodecon = cell(actnodes,3);

pores_no    = n/2*(m-1)+(n/2-1)*fix(m/2);                                  % n should be even (better to be factor of 4, m also even for this formula
pore_pos    = zeros(pores_no,4);

% ver 1.6, to define the type of pore in 3rd colomn
% to define if pores connected to BCs (BC1 or BC2 or ...) in 4th colomn
pore_iter   = 0;

realX       = [repmat([porelength/2 porelength],1,n/2) porelength/2];
node_realX  = zeros(actnodes,1);
node_realY  = zeros(actnodes,1);

stateoddX   = 1;

%%=========================================================================
%%-- Zentraler L"oser -----------------------------------------------------
%%=========================================================================

for i = 1:n                                                                % represents x axis
    if i/2 ~= round(i/2)                                                   % odd x
        for j = 1:m/2                                                      % represents y axis
            
            if stateoddX==1; 
                k          = 2*j-1;
                stateevenX = 2;
            else
                k          = 2*j;
                stateevenX = 4;
            end                                                            %k is actual y-axis position of active node
            actnodeno            = j + m/2*(i-1);
            nodes(k,i)           = actnodeno;
            x_actnode(actnodeno) = i;
            y_actnode(actnodeno) = k;
            
            %------------real position of nodes ---------------------------
            node_realX(actnodeno) = sum(realX(1:i-1));
            node_realY(actnodeno) = (k-1)*sqrt(3)*porelength/2;
            %--------------------------------------------------------------
            
            yx_actnodecon{actnodeno,1} = [k+1,i+1];
            yx_actnodecon{actnodeno,2} = [k-1,i+1];
            yx_actnodecon{actnodeno,3} = [k,i-1];
        end
        
    else                                                                   % even x
        for j = 1:m/2
            if stateevenX==2; 
                k = 2*j;stateoddX=3;
            else
                k = 2*j-1;stateoddX=1;
            end                                                            %k is actual y-axis position of active node
            
            actnodeno  = j+m/2*(i-1);
            nodes(k,i) = actnodeno;
            x_actnode(actnodeno) = i;
            y_actnode(actnodeno) = k;
            
            %------------real position of nodes ---------------------------
            node_realX(actnodeno) = sum(realX(1:i-1));
            node_realY(actnodeno) = (k-1)*sqrt(3)*porelength/2;
            %--------------------------------------------------------------
            
            i2 = actnodeno;
            % because all nodes around are Fard nodes, we use the k=2*j-1 to
            % find the actual number of it for
            
            yx_actnodecon{actnodeno,1}=[k-1,i-1];%((---1---))
            
            if y_actnode(actnodeno)~=1
                pore_iter             = pore_iter+1;
                pore_pos(pore_iter,2) = i2;                                  % pore ended node, even x
                if stateevenX == 2
                    pore_pos(pore_iter,1) = (yx_actnodecon{i2,1}(1)+1)/2+m/2*(yx_actnodecon{i2,1}(2)-1);%j=(k+1)/2 from k=2*j-1
                else
                    pore_pos(pore_iter,1) = yx_actnodecon{i2,1}(1)/2+m/2*(yx_actnodecon{i2,1}(2)-1);%j=k/2 from k=2*j
                end
                pore_pos(pore_iter,3) = 1;
            end
            
            yx_actnodecon{actnodeno,2}=[k+1,i-1]; %((---2---))
            
            if y_actnode(actnodeno)~=m
                pore_iter = pore_iter+1;
                pore_pos(pore_iter,2) = i2;                                % pore ended node, even x
                if stateevenX==2
                    pore_pos(pore_iter,1) = (yx_actnodecon{i2,2}(1)+1)/2 + m/2*(yx_actnodecon{i2,2}(2)-1);%j=(k+1)/2 from k=2*j-1
                else
                    pore_pos(pore_iter,1) = yx_actnodecon{i2,2}(1)/2 + m/2*(yx_actnodecon{i2,2}(2)-1);%j=k/2 from k=2*j
                end
                pore_pos(pore_iter,3) = 2;
            end
            
            yx_actnodecon{actnodeno,3} = [k,i+1];                          % ((---3---))
            
            if x_actnode(actnodeno)~=n
                pore_iter = pore_iter+1;
                pore_pos(pore_iter,2) = i2;                                % pore ended node, even x
                if stateevenX==2
                    pore_pos(pore_iter,1) = yx_actnodecon{i2,3}(1)/2 + m/2*(yx_actnodecon{i2,3}(2)-1); %j=k/2 from k=2*j
                else
                    pore_pos(pore_iter,1) = (yx_actnodecon{i2,3}(1)+1)/2 + m/2*(yx_actnodecon{i2,3}(2)-1); %j=(k+1)/2 from k=2*j-1
                end
                pore_pos(pore_iter,3) = 3;
            end
            
        end
        
    end
    
end

%%=========================================================================
% we can remove "nodes" variable, was for test to see if network works

node_no = actnodeno;

%maybe we can delete yx_actnodecon to free memory and use x_actnode and
%y_actnode and node_pores_nodes instead from now on
clear yx_actnodecon

% Obtaining the number(id) of pores connected to a specific node.
% V1.3, finding nodes and pores connected to a specific node.
node_pores_nodes = zeros(actnodes, 7);

% V1.6 added, the last colomn for specifying if the node is BC or not.
maxporepernod = 0;                                                         % number of allowed pore per node to inc the speed

%%=========================================================================

for inod=1:actnodes
    if (x_actnode(inod)==1) || (x_actnode(inod)==n)
        if (y_actnode(inod)==1) || (y_actnode(inod)==m)
            maxporepernod=1;
        else
            maxporepernod=2;
        end
    else
        if (y_actnode(inod)==1) || (y_actnode(inod)==m)
            maxporepernod=2;
        else
            maxporepernod=3;
        end
    end
    
    %=========== determining if node is BC1 or BC2 or none ======
    if x_actnode(inod)==1
        node_pores_nodes(inod,7)=1;
    elseif x_actnode(inod)==n
        node_pores_nodes(inod,7)=2;
    end
    %============================================================
    
    if rem(x_actnode(inod),2)~=0  %odd x, first column seek, fard kind
        seekcol=1;
        anothercol=2;
        type_node=1;
    else                          %even x, 2nd column seek, zoj kind
        seekcol=2;
        anothercol=1;
        type_node=2;%zoj
    end
    
    i2=0;
    for i=1:pores_no  %here we determine the node kind (fard or zoj) in order to determine the order of its neibours
        if i2==maxporepernod, break; end
        if pore_pos(i,seekcol)==inod
            i2=i2+1;
            node_con=pore_pos(i,anothercol);
            if type_node == 1
                if y_actnode(node_con)<y_actnode(inod)
                    position=1;
                elseif y_actnode(node_con)==y_actnode(inod)
                    position=2;
                elseif y_actnode(node_con)>y_actnode(inod)
                    position=3; 
                end
            else
                if y_actnode(node_con)<y_actnode(inod)
                    position=1;
                elseif y_actnode(node_con)==y_actnode(inod)
                    position=3;
                elseif y_actnode(node_con)>y_actnode(inod)
                    position=2;
                end
            end
            node_pores_nodes(inod,position)=i;
            %cols 1-3 for pores
            node_pores_nodes(inod,position+3)=pore_pos(i,anothercol);
            %cols 4-6 for nodes connected to those pores
            
            %====V1.6, if pores connected to BCs ===========
            if (x_actnode(inod)==n)||(x_actnode(node_con)==n)
                pore_pos(i,4)=2;
            end
            if (x_actnode(inod)==1)||(x_actnode(node_con)==1)
                pore_pos(i,4)=1;
            end
            %===============================================
        end
    end
end

%%=========================================================================
%%-- Postprocessing -------------------------------------------------------
%%=========================================================================

%%=========================================================================
%  Ver 1.8 to find the real places of nodes (X,Y) based on pore length unit

node_realXY  = [node_realX node_realY];
pore_diff_XY = (node_realXY(pore_pos(:,1),:) - node_realXY(pore_pos(:,2),:))/2;
pore_realX   = node_realXY(pore_pos(:,2),1) + pore_diff_XY(:,1);
pore_realY   = node_realXY(pore_pos(:,2),2) + pore_diff_XY(:,2);
pore_realXY  = [pore_realX pore_realY];

%%=========================================================================
% Finding the SCHEMATIC position of the porecenter for ploting the fluxes

xpore = zeros(pores_no,1);
ypore = zeros(pores_no,1);

logx = x_actnode(pore_pos(:,2)) - x_actnode(pore_pos(:,1));
logy = y_actnode(pore_pos(:,2)) - y_actnode(pore_pos(:,1));

xpore(logx==-1) = x_actnode(pore_pos((logx==-1),2)) + 0.5;
xpore(logx==1)  = x_actnode(pore_pos((logx==1),1))  + 0.5;
% xpore(logx==0)  = x_actnode(pore_pos((logx==0),1));

ypore(logy==-1) = y_actnode(pore_pos((logy==-1),2)) + 0.5;
ypore(logy==1)  = y_actnode(pore_pos((logy==1),1))  + 0.5;
ypore(logy==0)  = y_actnode(pore_pos((logy==0),1));

% for differenet state that hexagonals have a paralel side to X not Y

clear logx logy

%%-- Plotting -------------------------------------------------------------

%contour(XI, YI, u, 20);

% hold on;

figure;
% scatter(pore_realXY(:,1), pore_realXY(:,2), 50, '.');
scatter(node_realXY(:,1), node_realXY(:,2), 50, '*');
axis equal;

%%=========================================================================

% min(node_realX)*1000
% min(node_realY)*1000
% max(node_realX)*1000
% max(node_realY)*1000
% pores_no
% node_no

NodeData = struct('Nodes', node_realXY, 'NumberOfNodes', node_no, 'XNode', x_actnode, 'YNode', y_actnode, 'NodesOfPores', node_pores_nodes);
PoreData = struct('PoreInfo', pore_pos, 'PorePos', pore_realXY, 'NumberOfPores', pores_no);
GeometryData = struct('NodeData', NodeData, 'PoreData', PoreData, 'GeometryCoeffs', GeometryCoeffs);

end