function Main

%% Overview
% This Matlab file is the main file for a series of files:
% (i)   getGeometry.m, where the skeleton of the network is build
% (ii)  getRadii.m, where random values are assigned to each pore
% (iii) getFlow.m, where the fluid flow is solved for the network and
% (iv)  getTransport.m, where the transport of several species is
% computed.
%
% Further details on the used methods can be found in 'A reactive transport
% modeling approach to simulate biogeochemical processes in pore structures
% with pore-scale heterogeneities', which can be accessed by
% <http://dx.doi.org/10.1016/j.envsoft.2011.10.010>
%
% Copyright (C) 2013  Mehdi Gharasoo, Falk He"sse.
% The terms of the license agreement can be found in the README.txt.
 

%% Generating the geometry

%%  
% The geometrical parameters are defined as follows: 
%
% * n, which should be a factor of 4 and represents the X axis
% * m, which should be even and represents Y axis and
% * $$ l_{pore} $$ , which represents the pore length.
%

n               = 120;                                                      
m               = 36;                                                       
l_pore          = 0.001;                                                    

GeometryCoeffs  = struct('xIncr', n, 'yIncr', m, 'LengthOfPore', l_pore);

%%
% In the first step the raw geometry is generated, having a hexagonal grid
% and uniform pores of length $$ l_{pore} $$ .
%

GeometryData    = getGeometry(GeometryCoeffs);

%%
% This results in a uniform hexagonal grid.
%

%% Assigning the random pore radii 
%

%%
% The parameters of the random fields are:
%
% * $$ \lambda_x $$ , which is the correlation length in the x-direction 
% * $$ \lambda_y $$ , which is the correlation length in the y-direction
% * $$ \mu $$ , which is the expectation value of the random field and
% * $$ \sigma^2 $$ , which is the variance of the random field and
% * $$ func $$ , which defines the variogram function of the random field.
%

lambda_x        = 0.005;                                                   
lambda_y        = 0.005;                                                  
mu              = 0.16e-3;                                                 
sigma2          = 5e-9;
func            = 'Gauss';                                                    
RadiiCoeffs     = struct('lambda_x', lambda_x, 'lambda_y', lambda_y, 'mu', mu, 'sigma2', sigma2, 'func', func);                                 

%%
% In the second step random radii values $$ r_j $$ are assigned to 
% each pore $$ j $$ of the raw geometry.
%

RadiiData       = getRadii(GeometryData, RadiiCoeffs);             

%%
% After this step heterogeneous radii are assigned to each pore.
%

%% Calculating the flow field

%%
% The parameter of the flow field are
%
% * $$ \Delta p $$, which is the pressure drop along the geometry.

delta_p         = 62.15;
FlowCoeffs      = struct('delta_p', delta_p);

%%
% In the third step the volumetric flow $$ q_j $$ in each pore $$ j $$ is
% solved using the Hagen Poiseuille equation:
% $$ q_j = \frac{\pi r^4_j}{8 \nu} \frac{\Delta p_j}{l_{pore}} $$ . Here
% $$ \nu $$ is the dynamic viscosity of water and
% $$ \Delta p_j $$ is the pressure drop along pore $$ j $$ .
% Mass conservation is achieved by using Kirchhoff's law.
% The flow direction is assumed to be from the left to the right.
% The flow is used in order to determine the water velocity 
% $$ u_j = \frac{q_j}{\pi r_j^2} $$ in each pore $j$. These data are than
% used in order to determine the advective transport of each species.
%

FlowData        = getFlow(RadiiData, GeometryData, FlowCoeffs); 

%%
% After this calculation the flow values are assigned to each pore.
%

%% Computing the transport

%%
% The parameter of the transport are
%
% * $$ t_{end} $$, which is the lenght of the time interval,
% * $$ t_{0} $$, which is the beginning of the time interval,
% * $$ \Delta t $$, which is the lenght of the step size and
% * c_spec_no, which is the number of chemical species and
% * b_spec_no, which is the number of bacterial species.
%

t_end           = 30;
t_0             = 0;
dt              = 0.05;
c_spec_no       = 4;
b_spec_no       = 0;

TransportCoeffs = struct('t_end', t_end, 't_0', t_0, 'dt', dt, 'c_spec_no', c_spec_no, 'b_spec_no', b_spec_no);   

%% 
% In the last step the reactive transport for $$ N_{spec} $$ chemical
% species is solved using the advection-diffusion-reaction equation within
% each pore $$ j $$
% $$ \frac{\partial}{\partial t} c_{i,j} + \nabla \cdot ( u_j c_{i,j}) =
% D_i \Delta c_{i,j} + R_{i,j} $$ .
% Here $$ t $$ is the time, $$ c_{i,j} $$ is the concentration of species
% $$ i $$ in pore $$ j $$ , $$ D_i $$ is the diffusion cofficient of each
% species and $$ R_{i,j} $$ respective production and consumption rate.
% The last term is computed by virtue of the BRNS libary.
%

TransportData   = getTransport(FlowData, GeometryData, TransportCoeffs);
TransportDataV   = getTransportVanilla(FlowData, GeometryData, TransportCoeffs);


%%
% After this step concentration values for each species are assigned to each
% pore in each time step.
%

%% Postprocessing
% Several routines can be preformed here like plotting and saving.
%

poreXY = GeometryData.PoreData.PorePos;

% figure(2);
% scatter(poreXY(:,1), poreXY(:,2), 50, TransportData.value(:,end, 1), '.');
% axis equal;

data_diff = TransportData.value(:,end, 1) - TransportDataV.value(:,end, 1);
data_err = data_diff./TransportDataV.value(:,end, 1);

figure(2);
scatter(poreXY(:,1), poreXY(:,2), 50, TransportData.value(:,end, 1), '.');
scatter(poreXY(:,1), poreXY(:,2), 50, data_diff, '.');
axis equal;

% figure(2);
% scatter(poreXY(:,1), poreXY(:,2), 50, TransportDataV.value(:,end, 1), '.');
% axis equal;


end