% Calculates the growth rates for the osmotrophs. 
function mu = growthRates(mu,limFrac,limDepth)
%GROWTHRATES is a function calculates the growth rates in equations S16 to
%S18, in the Supporting Information.
%
%   Input:
%
%   mu - A structure containing storage vectors for the growth rates B, A,
%   and D, as well as the maximum growth rates for B, A, D, H, C, Z, and F.
%
%   limFrac - A structure containing the limitation fractions from the
%   equations S16 to S18. 
%
%   limDepth - A structure containing the logical nutrient and light
%   limitation depth indices. 
%
%   Output: 
%
%   mu - The structure containing the updated growth rates for B, A, and D,
%   as well as the constant maximum growth rates for B, A, D, H, C, Z, and
%   F.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% The bacterial growth rate: 
mu.B(limDepth.BP)  = mu.Bmax(limDepth.BP).*limFrac.BP(limDepth.BP); 
mu.B(limDepth.BL) = mu.Bmax(limDepth.BL).*limFrac.BL(limDepth.BL); 

% The autotrophic flagellate growth rate: 
mu.A(limDepth.AP)  = mu.Amax((limDepth.AP)).*limFrac.AP(limDepth.AP); 
mu.A(limDepth.AE) = mu.Amax(limDepth.AE).*limFrac.AE(limDepth.AE); 

% The diatom growth rate:
mu.D(limDepth.DP)  = mu.Dmax(limDepth.DP) .*limFrac.DP(limDepth.DP);
mu.D(limDepth.DS)  = mu.Dmax(limDepth.DS) .*limFrac.DS(limDepth.DS);
mu.D(limDepth.DE)  = mu.Dmax(limDepth.DE) .*limFrac.DE(limDepth.DE);

end
