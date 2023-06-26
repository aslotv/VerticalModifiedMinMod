% Performes an embedded explicit Runge-Kutta time step. 
function [u,e,growth] = eERK(un,dt,c,A,b,d,FSAL,sigma,mu,Y,I,H, ...
    limFrac,limDepth,f,delta,k,nu,K,P,R,rho,E_0,kappa,S_b, ...
    fishDistribution,FishIngestion,dz,zspan,N_z,N_vars)
% EERK is a function that performs an embedded explict Runge-Kutta step.
%
%   Input: 
% 
%   un - The solutions u at t = t_n. 
%   
%   dt - The time step.
% 
%   c - The nodes of the higher order method.
% 
%   A - The stage weights matrix.
% 
%   b - The scheme weigths of the lower order method.
% 
%   d - The difference between the scheme weights of the higher order
%   method and the scheme weights of the lower order method embedded
%   within.
%
%   FSAL - A logical, which equals true if the embedded Runge-Kutta scheme
%   has the "first same as last" property, false if not. 
% 
%   sigma - The mesozooplankton selectivity factor. 
% 
%   mu - A structure containing the maximum growth rates for B, A, D, H, C,
%   Z, and F, as well as storage vectors for the growth rates for B, A, and
%   D.
% 
%   Y - The yields. 
% 
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as storage vector for the ingestion rates, and the
%   denominator from the expression defining them.
% 
%   H - A structure containing the half-saturation constants for B, A, D,
%   H, C, Z, and F. 
% 
%   limFrac - A structure containing storage vectors for the limitation
%   fractions in equations S16 to S19 in the supplementary information. 
% 
%   limDepth - A structure containing storage vectors for logical index
%   vectors, where each non-zero elements corresponds to the indices where
%   the corresponding growth limitation is the most limiting. 
% 
%   f - The photosynthetic carbon overflow, the fraction of total loss that
%     enters detritus, and the fraction of total loss that is respired.
% 
%   delta - Specific mortality rates. 
% 
%   k - Dissolution, leakage and fragmentation rates, as well as empirical
%     coefficients related to chlorophyll light absorption.
% 
%   nu - A row vector of sinking rates, ordered according to the solution
%   matrices u and un.
% 
%   K - Attenuation from pure water and attenuation other than pure water
%   and chlorophyll.
% 
%   P - Photosynthetic quotient and inorganic phosphate concentration at
%   deep boundary (z_max).
% 
%   R - Respiratory quotient. 
% 
%   rho - The stoichiometric ratios and conversion factors.
% 
%   E_0 - The irradiance at the surface (z=0). 
% 
%   kappa - A row vector of turbulent diffusivities, ordered according to
%   the solution matrices u and un.
% 
%   S_b - The silicate concentration at the deep boundary (z = z_max)
% 
%   fishDistribution - The vertical distribution of fish.
% 
%   FishIngestion - The function handle for the function calculating the
%   fish ingestion rate stored in the variable expectedIngestion. 
% 
%   dz - The length of the depth steps.
% 
%   zspan - A column vector containing all the depths in the discretised
%   water column. 
% 
%   N_z - The number of depths in the discretised water column. 
% 
%   N_vars - The number of variables/equaitons in the system of PDEs. 
%
%   Output: 
% 
%   u - The approximated solutions.
% 
%   e - The proxy for the local error.
% 
%   growth - The growth in osmo- and phagotrophs.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

v = length(c); % Number of RK-stages.

% Storage array for the v increments in the time step:
Xi = zeros(N_z,N_vars,v);

% Storage matrices for the right hand side evaluations and the increments
% of the gross growth (i.e., growth without subtraction of loss):
dudt = zeros(N_z,N_vars,v);
vGrowth = zeros(N_z,7,v); 

% Calculating the v increments in the time step:
for i = 1:v
    
    % Calculating xi(i): 
    Xi(:,:,i) = un + dt*sum(A(i,:,1:i-1).*dudt(:,:,1:i-1),3);
   
    % Calculating the right hand side of the semi discretised PDEs as a
    % function of xi(i) and retrieving the growth:
    [dudt(:,:,i), vGrowth(:,:,i)] = RightHandSide(Xi(:,:,i),sigma,mu, ...
        Y,I,H,limFrac,limDepth,f,delta,k,nu,K,P,R,rho,E_0,kappa,S_b, ...
        fishDistribution,FishIngestion,dz,zspan,N_z,N_vars);
end

% Making a time step by defining u as the sum of un and dt times the linear
% combination of the time increments:
if FSAL
    % If the Runge-Kutta scheme has the First Same as Last (FSAL) property
    % then u is equal to xi_v:
    u = Xi(:,:,v); % In this case un + dt*sum(b.*dudt,3) = Xi(:,:,v). 
else
    % If the scheme does not have the FSAL property, then the linear
    % combination needs to be computed:
    u = un + dt*sum(b.*dudt,3);
end

% Calculating the proxy for the local error:
e = dt*sum(d.*dudt,3);

% The growth for the current time step is the linar combination of the
% growth in the v increments: 
growth = sum(b.*vGrowth,3);

end
