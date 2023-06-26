% Finds the logical limitation depth indices from the limitation fractions
% in the equations S16 to S18. 
function [limFrac,limDepth] = growthLimitation(un,E,H,limFrac,limDepth)
%LIMITATIONDEPTHINDICES This function calculates the limitation fractions
%in equations S16 to S18, in the Supporting Information, and finds the
%logical nutrient and light limitation depth indices. 
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   E - The irradiance for each depth in the discretised water column. 
%
%   H - The half-saturation constants for B, A, D, H, C, Z, and F.
%
%   limFrac - A structure with storage vectors for the limitation
%   fractions. 
%
%   limDepth - A structure with storage vectors for the limitation depth
%   indices. 
%
%   Output: 
%
%   limFrac - A structure with the limitation fractions. 
%
%   limDepth - A structure with the limitation depth indices.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes.  

% Calculating the limitation fractions in S16:
limFrac.BP = un(:,11)./(H.BP + un(:,11)); 
limFrac.BL = un(:,12)./(H.BL + un(:,12));

% The logical limitation depth indices for bacteria: 
limDepth.BP = limFrac.BP <= limFrac.BL;
limDepth.BL = ~limDepth.BP; 

% Calculating the limitation fractions in S17:
limFrac.AP = un(:,11)./(H.AP + un(:,11)); 
limFrac.AE = E./(H.AE + E);

% The logical limitation depth indices for autotrophic flagellates: 
limDepth.AP = limFrac.AP <= limFrac.AE;
limDepth.AE = ~limDepth.AP; 

% Calculating the limitation fractions in S18:
limFrac.DP = un(:,11)./(H.DP + un(:,11)); 
limFrac.DS = un(:,13)./(H.DS + un(:,13));
limFrac.DE = E./(H.DE + E);

% The logical limitation depth indices for diatoms: 
limDepth.DP = limFrac.DP <= min(limFrac.DS,limFrac.DE);
limDepth.DS = limFrac.DS < limFrac.DP & limFrac.DS <= limFrac.DE;
limDepth.DE = limFrac.DE < min(limFrac.DP,limFrac.DS); 

end
