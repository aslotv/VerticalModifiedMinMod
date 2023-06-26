function [dudt,growth,mu,limDepth,E,Chl,I,expectedIngestion] = ...
    RightHandSide(un,sigma,mu,Y,I,H,limFrac,limDepth,f,delta,k,nu,K, ...
    P,R,rho,E_0,kappa,S_b,fishDistribution,FishIngestion, ...
    dz,zspan,N_z,N_vars)
%RIGHTHANDSIDE This function calculates the right hand side of the partial
% differential equations S1-S15 in the Supporting Information, with
% the light regime, the vertical mixing, and the gravitational
% sinking accounted for.
%
%   Input: 
% 
%   un - The solutions u at t = t_n. 
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
%   dz - The length of the depth steps.
% 
%   zspan - A column vector containing all the depths in the discretised
%   water column. 
% 
%   N_z - The number of depths in the discretised water column. 
% 
%   N_vars - The number of variables/equaitons in the system of PDEs. 
% 
%   FishIngestion - The function handle for the function calculating the
%   fish ingestion rate stored in the variable expectedIngestion. 
%
%   Output:
%
%   dudt - The partial derivative of un with respect to time t.
%
%   growth - The growth in B, A, D, H, C, Z, and F.
% 
%   mu - A structure containing the maximum growth rates for B, A, D, H, C,
%   Z, and F, as well as the growth rates for B, A, and D.
% 
%   limDepth - A structure containing logical index vectors, where each
%   non-zero elements corresponds to the indices where the corresponding
%   growth limitation is the most limiting.
%
%   E - The irradiance for each depth in the discretised water column. 
%
%   Chl - Depth specific chlorophyll concentration. 
% 
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as the ingestion rates, and the denominator from the
%   expression defining them.
%
%   expectedIngestion - The expected value of the ingestion rate for fish
%   in the water column, according to the distribution fishDistribution, if
%   the fish are subject to diel vertical migration (DVM), or simply
%   returns the ingestion rate for fish if not.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Using the function ingestionRates to calculate the ingestion
% rates for H, C, Z and F:
I = ingestionRates(un,I,H,sigma);

% Calculates the fish ingesestion as the expected value of I_F
% according to the distribution fishDistribution if the fish are
% subject to diel vertical migration (DVM), or simply returns the
% ingestion rate I.F if not:
expectedIngestion = FishIngestion(zspan,I,fishDistribution);

% Using the irradiance function to calculate the water column
% irradiance, and retrieve the depth specific chlorophyll concentration:
[E,Chl] = irradiance(un,rho,K,k,E_0,zspan);

% Using the function growthLimitation to calculate the limitation fractions
% from equations S16 to S18, and retrieving the logical limitation indices: 
[limFrac,limDepth] = growthLimitation(un,E,H,limFrac,limDepth);

% Using the function growthRates to calculate the growth rates: 
mu = growthRates(mu,limFrac,limDepth);

% Using the function Growth to calculate the growth in bacteria, a.
% flagellates, diatoms, h. flagellates, ciliates, mesozooplankton,
% and fish:
growth = Growth(un,mu,Y,I,expectedIngestion,N_z);

% Using the function metabolicLossAllocation to calculate the
% allocation of metabolic losses:
[P,R,DOC_HCZ,DOC_F] = metabolicLossAllocation(un,growth, ...
    Y,I,expectedIngestion,P,R,f,rho);

% Using the function netGrowth to calculate the biological sources and
% sinks term, and to retrieve the selective ingestion proportionality
% factors:
s = netGrowth(un,growth,mu,Y,I,f,delta,k,P,R,rho,DOC_HCZ,DOC_F);

% Shifts the depth elements in u one step up, and one step down,
% while also incorporating the boundary conditions:
[u_shiftUpDiffusion,u_shiftDownDiffusion,u_shiftDownSinking] ...
    = shift(un,s,N_vars,P,S_b,kappa);

% Defining the right hand side:
dudt = s ...
    + dz.^(-2).*(kappa.*(u_shiftUpDiffusion - 2*un + u_shiftDownDiffusion))...
    - dz.^(-1).*(nu.*(un - u_shiftDownSinking));
end

