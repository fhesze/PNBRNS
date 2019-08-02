function FieldData = RandomField(GeometryData, coeffs, func)

% MATLAB script randomfield.m
%
% written by 
%
% Olaf. A. Cirpka
% Institute of Hydraulic Engineering
% University of Stuttgart
% Pfaffenwaldring 61
% 70550 Stuttgart
% Germany
% Olaf.Cirpka@IWS.Uni-Stuttgart.DE
%
% April 24, 2002 
%
% This script generates random periodic fields based on the
% parameters of a covariance function
%
% It requires:
% - the grid spacing of the unit cell
% - the type of covariance model 
% - the corresponding variance, correlation lengths and orientation
%
%
% E N J O Y   ! ! !
%

% clear all
% close all

%%=========================================================================
%%-- Defining Parameters --------------------------------------------------
%%=========================================================================

n     = GeometryData.GeometryCoeffs.xIncr;
m     = GeometryData.GeometryCoeffs.yIncr;
dPore = GeometryData.GeometryCoeffs.LengthOfPore;                          % L"ange einer Pore (in m)

X = 0.75.*dPore.*n;                                                        % Ausdehnung in x-Richtung (in m)
dx = 0.0001;                                                               % Schrittweite in x-Richtung (in m) 
XL = X/dx;                                                                 % totale Ausdehnung
x = (-XL/2:1:(XL-1)/2).*dx;                                                % 

Y = sin(1/3*pi)*dPore*m;                                                   % Ausdehnung in y-Richtung (in m)
dy = 0.0001;                                                               % Schrittweite in y-Richtung (in m) 
YL = Y/dy;                                                                 % totale Ausdehnung
y = (-YL/2:1:(YL-1)/2).*dy ;                                               % 

% geometrical parameters
[XI, YI] = meshgrid(x, y);                                                 % Grid in Physical Coordinates

% parameters of the field
mu      = coeffs(3);                                                       % Mittelwert
sigma2  = coeffs(4);                                                       % Standardabweichung

% misc
lx = [coeffs(1) coeffs(2)];
nx = [length(x) length(y)]; 
ntot = prod(nx);                                                           % total number of nodes 

%%=========================================================================
%%-- Field Generation Block -----------------------------------------------
%%=========================================================================

u = Fourier2D(XI, YI, [lx nx mu sigma2], func);

%%=========================================================================
%%-- Post Processing ------------------------------------------------------
%%=========================================================================

% to remove the minus cartezian position and put the origin on [0 0]
XImove = XI - XI(1,1);
YImove = YI - YI(1,1);

%%=========================================================================
%%-- Plotting -------------------------------------------------------------
%%=========================================================================

% scatter(XIvec,YIvec);

% contour(XI, YI, u);
pcolor(XImove, YImove, u);
shading flat;
daspect([1 1 1]);

%%=========================================================================
%%-- Saving ---------------------------------------------------------------
%%=========================================================================

% uncomment the following line to save the realization
% eval(sprintf('save field%3.3i.mat dx nx u -v4;', n));

FieldData = struct('value', u, 'xPos', XImove, 'yPos', YImove);

end

function u = Fourier2D(X, Y, params, func)

lx(1) = params(1);
lx(2) = params(2);
nx(1) = params(3);
nx(2) = params(4);

ntot = prod(nx);

mu = params(5);
sigma2 = params(6);

% ============== Auto-Covariance Block ====================================

H = sqrt((X/lx(1)).^2+(Y/lx(2)).^2);

% Covariance Matrix of Log-Conductivities
switch func
    case 'Exp'
        RYY = exp(-abs(H))*sigma2;
    case 'Gauss'
        RYY = exp(-H.^2)*sigma2;
    otherwise
        disp(['BEWARE:', func, ' is an unknown correlation function.']);
end

% ============== Power-Spectrum Block =====================================
% Fourier Transform (Origin Shifted to Node (1,1))

SYY = fftn(fftshift(RYY))/ntot;                                            % Yields Power Spectrum of the field
SYY = abs(SYY);                                                            % Remove Imaginary Artifacts
SYY(1,1,1) = 0;  
    
% nxhelp is nx with the first two entries switched
nxhelp = nx;

if(size(nx,2)>1)
    nxhelp(1:2)=[nx(2) nx(1)];
else
    nxhelp = [1,nx(1)];
end
    
% Generate a field of random real numbers,
% transform the field into the spectral domain,
% and evaluate the corresponding phase-spectrum.
% This random phase-spectrum and the given power spectrum
% define the Fourier transform of the random autocorrelated field.
u = sqrt(SYY).*exp(1i*angle(fftn(rand(nxhelp))));

% Backtransformation into the physical coordinates
u = real(ifftn(u*ntot)) + mu;

end
