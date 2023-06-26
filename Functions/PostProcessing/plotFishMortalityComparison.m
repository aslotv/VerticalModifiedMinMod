function [figure,tilingLayout,ax] = ...
    plotFishMortalityComparison(BaselineFileName,FishMortalityFileNames)
%PLOTFISHMORTALITYCOMPARISON This function makes a figure of tiled plots,
%each tile containing line plots of the vertical profiles of h. flagellates
%(H), ciliates (C), mesozooplankton (Z), bacteria (B), a. flagellates (A),
%diatoms (D), labile dissolved organic carbon (L), dissolved inorganic
%phosphate (P), and silicate (S), for different values of the
%FishMortalityFactor.
%
%   Input:
%
%   BaselineFileName (BaselineFileName = 'Baseline_steady.mat' as default)
%   - The name of the file containing the baseline parameters. Can also
%   explicitely request default value of BaselineFileName by giving
%   'default', 'Default', 'DEFAULT' (a case insensitive string comparison
%   is used), or the empty array [] as first input argument.
%
%   FishMortalityFileNames ('HalvedFishMortality_steady.mat' and
%   'TwofoldFishMortality_steady.mat' as default) - The names of the files
%   containing the simulations for the different values of specific
%   mortality rates for fish. This is a repeating argument type, give the
%   filenames as separate inputs (e.g.,
%   plotFishMortalityComparison(BaselineFileName,FishMortalityFileName1,...
%            FishMortalityFileName2,FishMortalityFileName3,...))
%
%   Output
%
%   This function has no output. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%% Arguments blocks
arguments
    BaselineFileName {mustBeFileOrDefaultOrEmpty(BaselineFileName)} ...
        = 'Baseline_steady.mat';
end

arguments(Repeating)
    FishMortalityFileNames {mustBeFileName}
end

%% Defining default values if nessecary

% Defining the default values of the names of the files containing the
% simulations with different values of the specific mortality rate for
% fish:
if isempty(FishMortalityFileNames)
    FishMortalityFileNames = { ...
        'HalvedFishMortality_steady.mat', ...
        'TwofoldFishMortality_steady.mat'...
        };
end

% Assigning the default value to BaslineFileName if BaselineFilename is
% equal to "default" or 'default' (regardless of the case of the letters)
% or if BaselineFileName = []:
if strcmpi(string(BaselineFileName),"default") || isempty(BaselineFilename)
    BaselineFileName ='Baseline_steady.mat';
end

%% Loading and storing solution structures and fish mortality factors

% Defining the baseline value for the fish mortality:
load(BaselineFileName,'Parameters')
baselineFishMortality = Parameters.delta.F;
clear Parameters

% Defining a storage vector for the fish mortalities:
FishMortalities = zeros(length(FishMortalityFileNames),1);

% Loading and storing fish mortalities and solution structures:
for i = 1:length(FishMortalityFileNames)
    temp = load(FishMortalityFileNames{i},'sol','Parameters');
    FishMortalities(i) = temp.Parameters.delta.F;
    sol(i) = temp.sol; %#ok<AGROW>
end

% Sorting the fish mortalities in ascending order and storing the
% sorting indices:
[FishMortalities,sortedIndices] = sort(FishMortalities);

% Using the sorting indices to sort the solution structures
% according to the sorted fish mortalities:
sol = sol(sortedIndices);

% Defining the fish mortality factors:
FishMortalityFactors = FishMortalities./baselineFishMortality;

%% Set up for the comparison plot

% The titles of the tiled plots:
Titles = {'H-Flag (nM-P)','Ciliates (nM-P)', 'Meso-Zooplankton (nM-P)', ...
    'Bacteria (nM-P)', 'A-Flag (nM-P)', 'Diatoms (nM-P)', ...
    'L-DOC (nM-C)','DIP (nM-P)', 'Silicate (nM-Si)'};

varIndices = [4:6 1:3 12 11 13]; % The indices of the variables in u.

% Storage cells:
PlotData = cell(1,length(varIndices));
x = PlotData;
Legend = cell(1,length(FishMortalityFactors));

% Defining the plot data and depth spans:
for j = 1:length(FishMortalityFactors)
    for i = 1:length(varIndices)
        PlotData{i}{j} = sol(j).u(:,varIndices(i));
        x{i}{j} = sol(j).zspan;
    end
    Legend{j} = [num2str(FishMortalityFactors(j)) '\delta_F'];
end

%% Making the figure
[figure,tilingLayout,ax] = plotInTiles(PlotData,x, ...
    "DepthPlot",true, ...
    'Layout',{3,3}, ...
    'FigureName','FishMortalityComparison', ...
    'Titles',Titles, ...
    'MainTitle','Vertical Profiles for Different Fish Mortality Factors', ...
    'MainXLabel','Concentrations (nM)', ...
    'MainYLabel','Depth (in meters)');

% Adding shared legend to the right side of the plots:
Legend = legend(Legend{:});
Legend.Layout.Tile = 'east';

%% Saving the figure

% Initial name of figure:
path = correctSlashes('SimulationResults\FinalResults\');
destination = 'FishMortalityComparison.fig';

answer = 'No'; % Initiates the while loop.
i = 0; % Initial value for the counter.

% As long as there already exists a file by the name in destination, a
% new name will be made before saving:
while isfile([char(path) destination]) && answer == "No"
    question = sprintf(['A figure-file by the name %s already exist, ' ...
        'do you wish to overwrite this file?'],destination);
    answer = questdlg(question,'File already exist','Yes','No','No');
    i = i+1; % Increases the counter by 1.
    switch answer
        case 'Yes'
           break
        case 'No'
            destination = ['(' num2str(i) ')FishMortalityComparison.fig'];
    end
end
savefig(figure,[path destination]);
end

function mustBeFileOrDefaultOrEmpty(BaselineFileName)
% Validation function for first input argument.

% True if input BaselineFileName is neither "default" nor "Default":
notDefault = BaselineFileName ~= "default" & BaselineFileName ~= "Default";

if ~isfile(BaselineFileName) && notDefault && ~isempty(BaselineFileName) ...
        && isempty(which(BaselineFileName))
    errID = 'MATLAB:validators:mustBeFile';
    errMSG = sprintf(['Invalid argument at position 1. ' ...
        'The following files do not exist: %s.'],BaselineFileName);
    throwAsCaller(MException(errID,errMSG));
end

end

function mustBeFileName(FileName)
% This function simply corrects the slashes (\ or /) in the file name given
% as input, if necessary, and passes the file name to the MATLAB(R)
% argument validation function mustBeFileName. If mustBeFile throws an
% error, the function checks if the file is in a folder added to the
% MATLAB(R) search path; if it is, then the argument is valid, if not, the
% error produced by mustBeFile is thrown. 
try
    FileName = correctSlashes(FileName);
    mustBeFile(FileName)
catch ME
    if isempty(which(FileName))
    throwAsCaller(ME);
    end
end
end
