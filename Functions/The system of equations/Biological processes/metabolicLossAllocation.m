% Calculates the allocation of the metabolic losses.
function [P,R,DOC_HCZ,DOC_F] = metabolicLossAllocation(un,growth, ...
    Y,I,expectedIngestion,P,R,f,rho)
%METABOLICLOSSALLOCATION is a function that calculates the metabolic losses
%from equations S23 to S33, in the Supporting Information.
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   growth - The growth in osmo- and phagotrophs. 
%
%   Y - The yields. 
%
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as the ingestion rates, and the denominator from the
%   expression defining them.
%
%   expectedIngestion - The expected value of the ingestion rate for F in the
%   water column. 
%
%   P - Photosynthetic quotient, inorganic phosphate concentration at deep
%   boundary (z_max).
%   
%   R - Respiratory quotient.                                                           
%
%   f - The photosynthetic carbon overflow, the fraction of total loss that
%   enters detritus, and the fraction of total loss that is respired.
%
%   rho - The stoichiometric ratios and conversion factors.
%
%
%   Output: 
%
%   P - Photosynthetic quotient, inorganic phosphate concentration at deep
%   boundary (z_max), and particulate losses.
%   
%   R - Respiratory quotient, and respiration losses.    
%
%   DOC_HCZ - DOC losses to L, from H, C, and Z. 
% 
%   DOC_F - DOC loss to L, from F. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%%% Total losses:
T_HC = (1-Y.H).*I.H.*un(:,4) + (1-Y.C).*I.C.*un(:,5); % From H and C.
T_Z  = (1-Y.Z).*I.Z.*un(:,6); % From Z.
T_F  = (1-Y.F).*expectedIngestion.*un(:,7); % From F.

%%% Particulate losses to detritus:
P.HC = f.d.*T_HC; % From H and C to Det_s.
P.Z  = f.d.*T_Z; % From Z to Det_f.
P.F = f.d.*T_F; % From F to Det_f.

%%% Respiration losses:
R.B   = (1/Y.BL-rho.B).*growth(:,1); % From B.
R.HCZ = f.rHCZ.*(T_HC + T_Z).*rho.CP; % From H, C, and Z.
R.F   = rho.CP.*f.rF.*T_F; % From F.

%%% DOC losses to L:
DOC_HCZ = (1 - f.rHCZ - f.d).*(T_HC + T_Z).*rho.CP; % From H, C, Z.
DOC_F   = rho.CP.*(1 - f.rF - f.d).*T_F; % From F.
end
