function NetworkData = MatchPoreNetwork(RandomField, GeometryData)

%% matching between generated pore network and a randm spatial distributer

poreXY  = GeometryData.PoreData.Pores;
poreNo  = GeometryData.PoreData.NumberOfPores;
radi    = zeros(poreNo,1);

% vectorizing the matrices
xVec    = reshape(RandomField.xPos,[],1);
yVec    = reshape(RandomField.yPos,[],1);
randVec = reshape(RandomField.value,[],1);                                 % random filed values

%%=========================================================================
%%-- Zentraler L"oser -----------------------------------------------------
%%=========================================================================

closest_p = ipdm(poreXY,[xVec yVec],'Subset','nearest','result','struct');

radi(closest_p.rowindex) = randVec(closest_p.columnindex);
radi_het = abs(radi);

%%=========================================================================
%%-- Plotting -------------------------------------------------------------
%%=========================================================================

% figure;
% scatter(poreXY(:,1), poreXY(:,2), 50, radi, '.');
% axis equal;
figure;
scatter(poreXY(:,1), poreXY(:,2), 50, radi_het, '.');
axis equal;

% max(poreXY(:,1))
% max(poreXY(:,2))
% min(poreXY(:,1))
% min(poreXY(:,2))
    
%save('..\radi_lx005_sigma2e-10.mat', 'radi_het')

NetworkData = struct('value', radi_het, 'xPos', xVec, 'yPos', yVec);

end