function expectedIngestion = FishIngestion(zspan,I,fishDistribution)
%FISHINGESTION This function estimates the expected value of the ingestion
%rate for fish in the discretised water column. 
%   
%   Input: 
%   
%   zspan - A vector containing all the discretized depths. 
%
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as the ingestion rates, and the denominator from the
%   expression defining them.
%
%   fishDistribution - The relative diel average distribution of fish in
%   the water column. 
%
%   Output: 
%   
%   expectedIngestion - The expected value of the ingestion rate for fish
%   in the discretised water column.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Calculates the *total ingestion by fish* in the water column:
expectedIngestion = trapz(zspan,I.F.*fishDistribution);

end

