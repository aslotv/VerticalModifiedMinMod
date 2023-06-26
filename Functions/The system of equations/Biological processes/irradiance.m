% Calculates the irradiance at each depth in zspan.
function [E,Chl] = irradiance(un,rho,K,k,E_0,zspan)
%IRRADIANCE is a function that calculates the irradiance in equation S35,
%in the Supporting Information.
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   rho - The stoichiometric ratios and conversion factors.
%
%   K - Attenuation from pure water and attenuation other than pure water
%   and chlorophyll.
% 
%   k - Dissolution, leakage and fragmentation rates, as well as empirical
%   coefficients related to chlorophyll light absorption.
%
%   E_0 - The surface irradiance. 
%
%   zspan - A row vector containing all the depths in the discretised water
%   column. 
%
%   Output: 
%
%   E - The irradiance for each depth in the discretised water column. 
%
%   Chl - Depth specific chlorophyll concentration. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Depth specific chlorophyll concentration:
Chl = rho.ChlP*(un(:,2)+un(:,3));

% Depth specific light attenuation:
K_z = cumtrapz(zspan,K.other + K.W + k.p*Chl.^k.e);

E = E_0*exp(-K_z); % Irradiance for each depth.
end
