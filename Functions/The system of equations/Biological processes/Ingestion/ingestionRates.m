% Calculates the ingestion rates.
function I = ingestionRates(un,I,H,sigma)
%INGESTIONRATES is a function that calculates the ingestion rates from
%equations S19 to S22, in the Supporting Information.
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   I - A structure containing storage vectors for the ingestion rates for
%   H, C, Z, and F, as well as the maximum ingestion rates for H, C, Z, and
%   F.
%
%   H - The half-saturation constants for B, A, D, H, C, Z, and F.
%
%   sigma - Mesozooplankton (Z) selectivity factor for C relative to D. 
% 
%   Output: 
%
%   I - A structure containing the ingestion rates for H, C, Z, and F, as
%   well as the maximum ingestion rates for H, C, Z, and F, and the
%   denominator from the expression defining them.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% The ingestion rate for heterotrophic flagellates: 
I.Hdenominator = H.H + un(:,1);
I.H = I.Hmax.*un(:,1)./I.Hdenominator;

% The ingestion rate for ciliates: 
I.Cdenominator = H.C + un(:,2) + un(:,4);
I.CA = I.Cmax.*un(:,2)./I.Cdenominator; % C ingestion of A. 
I.CH = I.Cmax.*un(:,4)./I.Cdenominator; % C ingestion of H. 
I.C = I.CA + I.CH; % Total C ingestion. 

% The ingestion rate for mesozooplankton: 
I.Zdenominator = H.Z + sigma*un(:,5) + un(:,3); 
I.ZD = I.Zmax.*un(:,3)./I.Zdenominator; % Z ingestion of D.
I.ZC = I.Zmax.*sigma.*un(:,5)./I.Zdenominator; % Z ingestion of C.
I.Z = I.ZD + I.ZC; % Total Z ingestion. 

% The ingestion rate for fish: 
I.Fdenominator = H.F + un(:,6);
I.F = I.Fmax.*un(:,6)./I.Fdenominator;

end
