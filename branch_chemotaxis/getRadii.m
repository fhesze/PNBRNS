function NetworkData = getRadii(GeometryData, coeffs)

%% Erzeugt Matrix der hydraulischen Leitf"ahigkeit K 

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

%%-- directory ------------------------------------------------------------

% home_path   = getSavePath;

%%-- deriving parameters --------------------------------------------------

nIncr       = GeometryData.GeometryCoeffs.xIncr;
mIncr       = GeometryData.GeometryCoeffs.yIncr;
dPore       = GeometryData.GeometryCoeffs.LengthOfPore;                    % L"ange einer Pore (in m)

poreXY      = GeometryData.PoreData.Pores;
nPore       = GeometryData.PoreData.NumberOfPores;

%%-- random field parameters ----------------------------------------------

lMaxX       = coeffs.lambda_x;                                             % maximum lenght in x direction 
lMaxY       = coeffs.lambda_y;                                             % maximum lenght in y direction 
lMinX       = 0;                                                           % minimum lenght in x direction 
lMinY       = 0;                                                           % minimum lenght in y direction 
H           = 0.33;                                                        % Hurst coefficient

muK         = coeffs.mu;                                                   % Expectation value of K
si2K        = coeffs.sigma2;                                               % Variance of K
muY         = log(muK) - 0.5.*log(1 + si2K./muK.^2);                       % Expectation value of Y
si2Y        = log(1 + si2K./muK.^2);                                       % Variance of Y
            
func        = coeffs.func;

%%-- numerical paramters ----------------------------------------------

nReal = 1;                                                                 % number of realizations
nMode = 2^(10);                                                            % number of modes

%%=========================================================================
%%-- Zentraler L"oser -----------------------------------------------------
%%=========================================================================

% u = zeros(Nx, Ny);

for n = 1:nReal
    
    Ksi = randn(nMode, 2);                                                 % gaussion-distributed random variable for sampling
    Nu  = getRandomSet([lMaxX lMaxY lMinX lMinY H nMode], func);         % user-specified random variable
  
    for i = 1:nPore
        u(i) = calculate_Y(poreXY(i,1), poreXY(i,2), nMode, Ksi, Nu);
    end
    
    Y   = sqrt(si2Y).*u + (muY);                                           % computing the log-hydraulic conductivity Y = ln(K) 
    K   = exp(Y);                                                          % computing the hydraulic conductivity K

end

%%=========================================================================
%%-- Postprocessing -------------------------------------------------------
%%=========================================================================

load('radi_lx005_sigma5e-9_50x50mm','radi_het');
radi        = radi_het;

%%-- Plotting -------------------------------------------------------------

%contour(XI, YI, u, 20);

figure;
scatter(poreXY(:,1), poreXY(:,2), 50, radi, '.');
axis equal;

%%-- Assembling data ------------------------------------------------------

NetworkData = struct('value', K', 'xPos', poreXY(:,1), 'yPos', poreXY(:,2));

end

function u = calculate_Y(XI, YI, n0, Ksi, Nu)

u = zeros(size(XI, 1), size(XI,2));

for j = 1:n0
    
    theta = (Nu(j,1).*XI + Nu(j,2).*YI);

    A = Ksi(j,1)*cos(theta);
    B = Ksi(j,2)*sin(theta);
    u = u + (1/sqrt(n0))*(A + B);
        
end
    
end

function Nu = getRandomSet(params, flag)

lMaxX       = params(1);
lMaxY       = params(2);
lMinX       = params(3);
lMinY       = params(4);
H           = params(5);
n0          = params(6);  

%%-- Erzeugung von Zufallssequenz Nu --------------------------------------

% uniformly-distributed random variable for the correlation function
gam1        = rand(n0, 1);
gam2        = rand(n0, 1);

Nu          = zeros(n0, 2);
    
switch flag
    case 'TruncGauss'
        r = getInverse2D([1 1 H], (2.*gam1 - 1), flag);
        Nu(:,1) = r.*cos(2.*pi.*gam2')./lMaxX;
        Nu(:,2) = r.*sin(2.*pi.*gam2')./lMaxY;
    case 'TruncExp'
        r = getInverse2D([1 1 H], (2.*gam1 - 1), flag);
        Nu(:,1) = r.*cos(2.*pi.*gam2')./lMaxX.*pi;
        Nu(:,2) = r.*sin(2.*pi.*gam2')./lMaxY.*pi;
    case 'Gauss'
        % Nu(:,1) = getInverse([corrLenX 1 H], (2.*gam1 - 1), flag)./lMaxX;
        % Nu(:,2) = getInverse([corrLenY 1 H], (2.*gam2 - 1), flag)./lMaxY;
%         r=sqrt(-log(1-gam1));
%         Nu(:,1) = lMaxX.*r.*cos(2.*pi.*gam2)./lMaxX;
%         Nu(:,2) = lMaxY.*r.*sin(2.*pi.*gam2)./lMaxY;
        Nu(:,1) = sqrt(2).*randn(n0, 1)./lMaxX;
        Nu(:,2) = sqrt(2).*randn(n0, 1)./lMaxY;
        % r = getInverse2D([corrLenX./sqrt(8/pi) 1 H], (2.*gam1 - 1), flag);
        % Nu(:,1) = r.*cos(2.*pi.*gam2');
        % Nu(:,2) = r.*sin(2.*pi.*gam2');
    case 'Exp'
        Nu(:,1) = sqrt(1./gam2.^2 - 1).*cos(2.*pi.*gam1)./lMaxX;
        Nu(:,2) = sqrt(1./gam2.^2 - 1).*sin(2.*pi.*gam1)./lMaxY;
        % nu1 = getInverse([corrLenX 1 H], (2.*gam1 - 1), flag);
        % nu2 = getInverse([corrLenY 1 H], (2.*gam2 - 1), flag);
end


    
end

function y = getInverse2D(params, x0, type)

    for i = 1:length(x0)
        
        y(i) = getInvScal(params, x0(i), type);
        
    end

end

function y = getInvScal(params, x0, type)

    lambda = params(1);
    C = params(2);
    H = params(3);
    H2 = 2.*H;
    
    switch type
        case 'Gauss' 
            y = fzero(@Gauss, [0.5]);
        case 'Exp' 
            % y = fzero(@Exp, [tan(x0.*pi./2)]);
        case 'TruncGauss'
            y = fzero(@TruncGauss2D, [tan(x0.*pi./2)]);
        case 'TruncExp'
            y = fzero(@TruncExp2D, [tan(x0.*pi./2)]);   
    end

    function f = Gauss(x)
        
        h = CumDistFunc2D(x, [lambda H], 'Gauss');
        f = h - x0;
        
    end
    
    function f = TruncGauss2D(x)
        
        h = CumDistFunc2D(x, [lambda H], 'TruncGauss');
        f = h - x0;
        
    end

    function f = TruncExp2D(x)
        
        h = CumDistFunc2D(x, [lambda H], 'TruncExp');
        f = h - x0;
        
    end

end

function CDF = CumDistFunc2D(k, params, flag)

lambda = params(1);
% lambda = lambda./pi;
H = params(2);

nUp = 1./params(1);
H2 = 2.*H;

switch flag
    case 'Gauss'
        % CDF = erf(k./(2).*lambda);
        lambda = lambda.*pi./2;
        CDF = sign(k).*(1 - exp(-(lambda*k/2).^2));
    case 'Exp'
        % not yet implemented        
    case 'TruncGauss'
        var = (k./nUp).^2./pi;
        A = exp(-var);
        B = var.^(-H).*(gamma(1 + H) - myGammaInc(1 + H, var));
        CDF = sign(k).*(1 - A - B);
    case 'TruncExp'
        var =  (pi.*k).^2;
        A = H./2.*gamma(1+H)./gamma(2+H).*(pi.*k).^2;
%         B = hypergeom([1.5, (1+H), 1], [(2+H), 2], -(pi.*k).^2);
%         B = hypergeomLaplace([1.5, (1+H), 1], [(2+H), 2], -(pi.*k).^2);
        B = ((2+H2).*(sqrt(1 + var) - 1)./var - hypergeomLaplace(1, 0.5+H, 2+H, -var))./(H.*sqrt(1 + var));
        CDF = sign(k).*A.*B;
end

CDF(find(k == 0)) = 0;

end

function gamma_inc = myGammaInc(a, x)

%% Definition of the incomplete Gamma function according to Mathematica

gamma_inc = gamma(a).*gammainc(x, a, 'upper');

end
