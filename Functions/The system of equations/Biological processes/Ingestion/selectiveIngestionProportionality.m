function  [selectProp_CA,selectProp_CH,selectProp_ZD, selectProp_ZC] ...
    = selectiveIngestionProportionality(un,sigma)
%SELECTIVEINGESTIONPROPORTIONALITY This function calculates the fractions
%representing loss by ingestion of predator, appearing in the equations S2
%to S5, in the Supporting Information, are defined in this function.
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   sigma - The mesozooplankton selectivity factor. 
%
%   Output: 
%
%   selectProp_CA - The fraction of A ingested, by C, of the total C
%   ingestion.
%
%   selectProp_CH - The fraction of H ingested, by C, of the total C
%   ingestion.
%
%   selectProp_ZD - The fraction of D ingested, by Z, of the total Z
%   ingestion.
%
%   selectProp_ZC - The fraction of C ingested, by Z, of the total Z
%   ingestion.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Proportionality constants of selective ingestion by ciliates (C)
% of autotrophic- and heterotrophic flagellates, (A) and (H)
% respectively:
selectProp_CA = un(:,2)./(un(:,2) + un(:,4)); % From S2
selectProp_CH = un(:,4)./(un(:,2) + un(:,4)); % From S4

% Proportionality constants of selective ingestion by
% mesozooplankton (Z) of diatoms (D) and ciliates (C):
selectProp_ZD = un(:,3) ...
    ./(un(:,3) + sigma*un(:,5)); % From S3
selectProp_ZC = (sigma*un(:,5)) ...
    ./(un(:,3) + sigma*un(:,5)); % From S5
end
