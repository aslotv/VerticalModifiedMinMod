function expectedIngestion = FishIngestionNoDVM(~,I,~)
%FISHINGESTIONNoDVM A function that simply returns a copy of the ingestion rate.
%   
%   Input: 
%
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as the ingestion rates, and the denominator from the
%   expression defining them.
%
%   Output: 
%   
%   expectedIngestion - The ingestion rate for fish.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Simply returns a copy of the ingestion rate I_F:
expectedIngestion = I.F;

end

