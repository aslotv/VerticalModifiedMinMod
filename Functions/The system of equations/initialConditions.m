% Defines the initial conditions.
function [u_t0, u0, growth0, u0Mean, TotalP0mean, Production0] ...
    = initialConditions(mode,options)
%INITIALCONDITIONS is a function that defines the initial conditions for
%the Vertically Modified MinMod.
%
%   Input:
%
%   mode - Must be either 'load' or 'set'. With mode='set' the initial
%   conditions are set to pre-defined values, with mode='load', the initial
%   conditions are set to the solutions from the last time step, loaded
%   from the file provided as optional input or from the default file if
%   no file name is provided.
%
%   Optinal input:
%
%   N_z (N_z = 701 as default) - Number of depths in zspan, given if
%   mode='set'.
%
%   N_vars (N_vars = 15 as default) - Number of state variables or
%   equations, given if mode='set'.
%
%   zspan - A vector of all the depths in the discretised water column. 
%
%   zrange - The depth range of the water column. 
%
%   FileName (FileName = 'Baseline_steady.mat' as default) - The name of
%   the file in which the results from the previous run is stored, given if
%   mode='load'.
%
%   ResetOxygen (ResetOxygen = false as default) - Will reset oxygen to the
%   initial value from t=0 if true, given if mode='load'.
%
%   NoFish (NoFish = false as default) - NoFish = true results in a
%   simulation in which fish is not included.
%
%   OxyObsFileName (OxyObsFileName = 'oxygenProfile_0to700m.mat' as
%   default) - The name of the MAT-file in which the oxygen observations
%   are stored.
%
%   OxyObsVarName (OxyObsVarName = 'oxy' as default) - The name of the
%   variable in the MAT-file containing the oxygen observations.
%
%   FishDistFileName (FishDistFileName = 'fishDistribution_0to700m.mat' as
%   default) - The name of the MAT-file in which the fish distribution is
%   stored.
%
%   FishDistVarName (FishDistVarName = 'fishDistribution' as default) - The
%   name of the variable containing the fish distribution.
%
%   ContinuedSimulation (ContinuedSimulation = true as default) - A logical
%   telling this function whether or not the current simulation is a
%   continuation of a simulation that produced the results in the loaded
%   MAT-file. 
%
%   Output
%
%   u_t0 - The initial conditions for the current simulation.
%
%   u0 - The initial conditions from the first simulation if mode = 'load'
%   or the solutions structure loaded does not have u0 as a field. If mode
%   = 'load' and u0 is a field in the solutions structure, then u0 is set
%   equal to the u0 from the loaded solutions structure (i.e., u0 =
%   sol.u0).
%
%   growth0 - The initial value for the osmo- and phagotroph growth. 
%
%   u0Mean - The initial average water column concentration. 
%
%   TotalP0mean - The initial average water column concentration of total
%   P.
%
%   Production0 - The initial values for the osmo- and phagrotrophs
%   production. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

arguments
    mode {mustBeMember(mode,{'load','set'})}
    options.N_z (1,1) {mustBeInteger,mustBePositive}
    options.N_vars (1,1) {mustBeInteger,mustBePositive}
    options.zspan (:,1) {mustBeVector}
    options.zrange (1,1) double {mustBeNonzero}
    options.FileName {mustBeFileName(options.FileName)}
    options.ResetOxygen (1,1) {mustBeNumericOrLogical}
    options.NoFish (1,1) {mustBeNumericOrLogical}
    options.OxyObsFileName {mustBeFileName(options.OxyObsFileName)}
    options.OxyObsVarName {mustBeTextScalar}
    options.FishDistFileName {mustBeFileName(options.FishDistFileName)}
    options.FishDistVarName {mustBeTextScalar}
    options.ContinuedSimulation {mustBeNumericOrLogical}
end

% Assigns a cell array containing the input given to this function in the
% workspace of the caller of this function:
initialConditionsInput = ...
    [{'mode',mode} ...
    reshape([fieldnames(options)'; struct2cell(options)'],1,[])];
assignin("caller","initialConditionsInput",initialConditionsInput);

% Defines default values for missing fields in options:
options = setDefaultOptions(options);

switch mode

    case 'set'

        % Initial conditions at surface:
        u_t0 = ones(1,options.N_vars);
        u_t0([11 13]) = [1000 1000*10.7]; % For P and S.

        % Initial conditions throughout water column:
        u_t0 = repmat(u_t0,options.N_z,1); % Homogeneous initial conditions.

        if options.NoFish
            % Sets initial conditions for fish to zero if simulation
            % without fish is desired:
            u_t0(:,7)=0;
        else
            % Fish does not have homogeneus initial conditions:
            fishBiomass = 3 ... from g WW m^(-2)
                *1e9*0.1*0.5/(12*106); % to nmol P m^(-2).

            % Loading vertical fish distribution:
            fishdist = ...
                load(options.FishDistFileName,options.FishDistVarName);
            fishDistribution = fishdist.(options.FishDistVarName);

            % Allocating initial biomass vertically.
            u_t0(:,7) = fishBiomass*fishDistribution/1000;
        end

        % Oxygen does not have homogeneus initial conditions:

        % Loading the oxygen observations and storing it in a
        % structure:
        OxyObs = load(options.OxyObsFileName,options.OxyObsVarName);

        % Extracting the variable with the oxygen observations from the
        % structure:
        oxy = OxyObs.(options.OxyObsVarName);

        % Initial values for oxygen is set equal to observations:
        u_t0(:,15) = oxy;

        clear oxy % Cleans up workspace.

        % Storing these initial conditions for keeping track for future
        % subsequent runs:
        u0 = u_t0;

        % The initial growth of the osmo- and phagotrophs: 
        growth0 = zeros(options.N_z,7); 

        % The initial average water column concentration and total P
        % concentration:
        u0Mean = trapz(options.zspan,u0)/options.zrange;
        TotalP0mean = sum(trapz(options.zspan,u0(:,1:11))) ...
            /options.zrange;

        % The initial production of the osmo- and phagotrophs: 
        Production0 = zeros(1,7); 

    case 'load'

        % The name of the MAT-file to load, including the path from the
        % current directory to the MAT-file with path separators corrected
        % if neccessary:
        MATfile  = correctSlashes(options.FileName);

        % Loading the contents of the MAT-file and storing it in a
        % structure:
        loadedVariables = load(MATfile);

        % Asserting that the MAT-file contains a solution structure (sol)
        % with a field u (the initial conditions for the current
        % simulation):
        assert(isfield(loadedVariables,'sol') && ...
            isfield(loadedVariables.sol,'u'),['The initial conditions ' ...
            'must be stored in a field called "u" in a solution ' ...
            'structure called "sol".']);

        % The variable u is extracted from the solution structure sol,
        % which is stored in the structe loadedVariables: 
        u_t0 = loadedVariables.sol.u; % Solutions from previous run.

        % Checking if the solutions structure loaded from the MAT-file has
        % the required fields: 
        [hasRequiredFields,missingFields] ...
            = checkFields(loadedVariables.sol); 

        % If the loaded solutions structure does not have the required
        % fields and the optional input argument ContinuedSimulation equals
        % true, then the current simulation is treated as a new simulation:
        if ~hasRequiredFields && options.ContinuedSimulation
            
            % The value of ContiuedSimulation is changed to false: 
            options.ContinuedSimulation = false; 
            
            % The cell array initialConditionsInput is updated:
            TF = cellfun(@(x) string(x)=="ContinuedSimulation", ...
                initialConditionsInput);

            if any(TF)
                % If ContinuedSimulation was given as an input argument,
                % then the value is updated within the cell array
                % initialConditionsInput: 
                 initialConditionsInput{find(TF==1)+1} = false; 
            else
                % If the default value was used, the updated value is added
                % to the cell array initialConditionsInput: 
                initialConditionsInput = [initialConditionsInput ...
                    {"ContinuedSimulation",false}]; 
            end
            
            % The updated cell array initialConditionsInput is reassigned
            % to the callers workspace: 
            assignin("caller","initialConditionsInput", ...
                initialConditionsInput);

            % Issuing a warning: 
            warning(['\nThe loaded solution structure does not have ' ...
                'all the required fields and the current simulation ' ...
                'will therefore be treated as a new simulation.' ...
                '\nThe loaded solution structure is missing the ' ...
                'following required fields: %s.\n'],missingFields)
        end

        % Initial conditions from t=0:
        if options.ContinuedSimulation % If continued simulation. 
            u0 = loadedVariables.sol.u0;

            % The initial growth of the osmo- and phagotrophs:
            growth0 = loadedVariables.sol.growth;

            % The initial average water column concentration and total P
            % concentration:
            u0Mean = loadedVariables.sol.uMean(end,:); 
            TotalP0mean = loadedVariables.sol.TotalPmean(end,:); 

            % The initial production of the osmo- and phagotrophs:
            Production0 = loadedVariables.sol.Production(end,:);

        else % If new simulation.
            u0 = u_t0;

            % The initial growth of the osmo- and phagotrophs:
            growth0 = zeros(options.N_z,7);

            % The initial average water column concentration and total P
            % concentration:
            u0Mean = trapz(options.zspan,u0)/options.zrange;
            TotalP0mean = sum(trapz(options.zspan,u0(:,1:11))) ...
                /options.zrange;

            % The initial production of the osmo- and phagotrophs:
            Production0 = zeros(1,7);
        end

        % If desired, reset values for oxygen to observed values:
        if options.ResetOxygen

            % Loading the oxygen observations and storing it in a
            % structure:
            OxyObs = load(options.OxyObsFileName,options.OxyObsVarName);

            % Extracting the variable with the oxygen observations from the
            % structure:
            oxy = OxyObs.(options.OxyObsVarName);

            % Initial values for oxygen is set equal to observations:
            u_t0(:,15) = oxy;

            clear oxy OxyObs % Cleans up workspace.
        end

        % If simulation without fish is desired, set initial conditions for
        % fish to zero:
        if options.NoFish; u_t0(:,7)=0; end

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

function options = setDefaultOptions(options)
% This function finds the options that wasn't specified as
% name-value argument pairs, and sets them to the default value.

% The names of all the fields in the options structure:
fieldNames = split("N_z, N_vars, zspan, zrange, FileName, " + ...
    "ResetOxygen, NoFish, OxyObsFileName, OxyObsVarName, " + ...
    "FishDistFileName, FishDistVarName, ContinuedSimulation",", ");

default.N_z = 701; % Default number of discretized depths.
default.N_vars = 15; % Default number of equations in the system.

default.zspan = 0:700; % The default depth span. 
default.zrange = 700; % The default depth range. 

% Default name of the file to be loaded if mode = 'load':
default.FileName = 'Baseline_steady.mat';

% Default value for the option to reset the values of the initial
% conditions for oxygen if mode = 'load':
default.ResetOxygen = false;

% Default value for the options to run the simulation with no fish:
default.NoFish = false;

% Default name of the file and variable containing the oxygen
% observations:
default.OxyObsFileName = 'oxygenProfile_0to700m.mat';
default.OxyObsVarName = 'oxy';

% Default name of the file and variable containing the fish
% distribution:
default.FishDistFileName = 'fishDistribution_0to700m.mat';
default.FishDistVarName = 'fishDistribution';

% Default value for the logical determining whether or not the current
% simulation is a continuation of a simulation that produced the results
% in the loaded MAT-file. 
default.ContinuedSimulation = true; 

% Assigns default values for the missing fields in the options
% structure:
for i = 1:length(fieldNames)
    if ~isfield(options,fieldNames(i))
        options.(fieldNames(i)) = default.(fieldNames(i));
    end
end
end


function [hasRequiredFields,missingFields] = checkFields(sol)

% The field names of the solution structure: 
fieldNames = fieldnames(sol);

% The required field names: 
requiredFieldNames = {'zspan','tspan','u0','uMean', ...
    'TotalPmean','Production','u','un','growth'}';

% A logical vector, whose non-zero elements correspond to the required
% field names that are in the solution structure (sol): 
TF = cellfun(@(x) any(strcmp(fieldNames,x)),requiredFieldNames); 

if all(TF)
    hasRequiredFields = true; 
else
    hasRequiredFields = false; 
end

% The names of the missing fields: 
missingFields = string(strjoin(requiredFieldNames(~TF))); 
end
