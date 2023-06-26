function [NoFish,NoDVM,DiffusivityFactor,FishMortalityFactor,Scenario] = ...
    SetScenario(Scenario)
%SETSCENARIO is a function that assigns values to the logical variables
%NoFish and NoDVM, and to the factors multiplied with the parameters kappa
%(turbulent diffusivity) and delta.F (specific mortality rate for fish),
%based on the chosen scenario. It is also possible to set these values
%manually.
%
%   Input:
%
%   Scenario - Must equal one of the following:
%
%            * 'Baseline' - Assigns values to output corresponding to the
%            baseline simulation, i.e., NoFish and NoDVM equals false, and
%            DiffusivityFactor and FishMortalityFactor equals 1.
%
%            * 'HalvedDiffusivity' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            DiffusivityFactor, which equals 0.5.
%
%            * 'TwofoldDiffusivity' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            DiffusivityFactor, which equals 2.
%
%            * 'TenfoldDiffusivity' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            DiffusivityFactor, which equals 10.
%
%            * 'HalvedFishMortality' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            FishMortalityFactor, which equals 0.5.
%
%            * 'TwofoldFishMortality' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            FishMortalityFactor, which equals 2.
%
%            * 'NoFish' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            NoFish and NoDVM, which both equals true.
%
%            * 'NoDVM' - Assigns values to the output
%            corresponding to the baseline simulation, except for
%            NoDVM, which equals true.
%
%            * 'SetManually' - Assigns values to the output
%            corresponding to the input given when prompted.
%
%   Output:
%
%   NoFish - NoFish = true results in a simulation without fish and without
%   DVM. NoFish = false results in a simulation with fish, however NoDVM
%   can be either true or false.
%
%   NoDVM - NoDVM = true results in a simulations where fish is
%   subjected to diel vertical migration; NoDVM = false results in a
%   simulation where fish is not subject to diel vertical migration.
%
%   DiffusivityFactor - The value set to DiffusivityFactor will be
%   multiplied with the baseline value of the turbulent diffusivity (kappa),
%   defined in the function "defineParameters.m".
%
%   FishMortalityFactor - The value set to FishMortalityFactor will be
%   multiplied with the baseline value of the specific mortality rate for
%   fish (delta.F), defined in the function "defineParameters.m".
%
%   Scenario - The text scalar Scenario, given as input.
%
%   Output usage:
%
%   NoFish - Used in functions "initialConditions.m" and
%   "defineParameters.m".
%
%   NoDVM - Used in functions "initialConditions.m" and
%   "saveAndPlotResults.m".
%
%   DiffusivityFactor - Used in the function "defineParameters.m".
%
%   FishMortalityFactor - Used in the function "defineParameters.m".
%
%   Scenario - Used in the function "saveAndPlotResults.m"
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

arguments
    Scenario {mustBeMember(Scenario,{'Baseline','HalvedDiffusivity', ...
        'TwofoldDiffusivity','TenfoldDiffusivity', ...
        'HalvedFishMortality','TwofoldFishMortality','NoFish','NoDVM', ...
        'SetManually'})}
end

% Defining the output arguments (save for Scenario) as empty arrays, in
% case function is canceled:
[NoFish,NoDVM,DiffusivityFactor,FishMortalityFactor] = deal([]);

switch Scenario
    case 'Baseline'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 1;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'HalvedDiffusivity'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 0.5;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'TwofoldDiffusivity'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 2;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'TenfoldDiffusivity'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 10;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'HalvedFishMortality'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 1;

        % Using baseline fish mortality:
        FishMortalityFactor = 0.5;

    case 'TwofoldFishMortality'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = false;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 1;

        % Using baseline fish mortality:
        FishMortalityFactor = 2;

    case 'NoFish'

        % F is included in the simulation:
        NoFish = true;

        % F is subject to diel vertical migration:
        NoDVM = true;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 1;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'NoDVM'

        % F is included in the simulation:
        NoFish = false;

        % F is subject to diel vertical migration:
        NoDVM = true;

        % Using baseline turbulent diffusivity:
        DiffusivityFactor = 1;

        % Using baseline fish mortality:
        FishMortalityFactor = 1;

    case 'SetManually'

        answer = questdlg('Include fish in the simulation?', ...
            'Inclution of fish in simulation',...
            'Yes','No','Yes');

        if isempty(answer)
            % If cancel button was pressed, end function: 
            return
        elseif answer == "Yes"

            % If user answers yes to include fish in the simulation, set
            % NoFish equal to false and ask if they want to run a
            % simulation with DVM:
            NoFish = false;
            answer = questdlg(['Subject fish to diel vertical ' ...
                'migration?'],'Diel Vertical Migration','Yes','No','Yes');

            % If user answers yes to run simulation with DVM, NoDVM
            % equals to false, if user answers no, NoDVM equals true:
            if ~isempty(answer)
                NoDVM = answer == "No";
            else
                % If cancel button was pressed, end function:
                return
            end

        elseif answer == "No"

            % If user answers no to include fish in the simulation, set
            % NoFish equal to true:
            NoFish = true;

            % In a simulation with no fish, the value of NoDVM is arbitrary
            % with respect to the simulation results. However, to avoid
            % uneccessary computaion, NoDVM is set equal to true:
            NoDVM = true;
        end

        wrongInput = true; % Equals true until user gives correct input.
        while wrongInput % While the user gives wrong input

            % Ask user to give input:
            DiffusivityFactor = ...
                inputdlg(['Enter desired value of DiffusivityFactor ' ...
                '(must be scalar) or press enter for default value ' ...
                'DiffusivityFactor'],'Diffusivity factor',1,{'1'});

            % Checks if the user clicked the cancel button:
            if isempty(DiffusivityFactor)
                % If cancel button was pressed, end function:
                return
            end

            % Checks if correct input was given for the diffusivity factor
            % and returns the input given as a double of length 1 if
            % correct:
            [wrongInput,DiffusivityFactor] = checkInput(DiffusivityFactor);
        end


        wrongInput = true; % Equals true until user gives correct input.
        while wrongInput
            FishMortalityFactor = inputdlg(['Enter desired value of ' ...
                'FishMortalityFactor (must be scalar) or press enter ' ...
                'for default value FishMortalityFactor = 1.'], ...
                'Diffusivity factor',1,{'1'});

            % Checks if the user clicked the cancel button:
            if isempty(FishMortalityFactor)
                return
            end

            % Checks if correct input was given for the diffusivity factor
            % and returns the input given as a double of length 1 if
            % correct:
            [wrongInput,FishMortalityFactor] ...
                = checkInput(FishMortalityFactor);
        end

        Scenario = string(inputdlg(['Enter an alternative prefix ' ...
            'to the steady state MAT-file (will be saved as ' ...
            'prefix_steady.mat) or press enter to use SetManuall as ' ...
            'prefix (SetManually_steady.mat)'],['Prefix to steady ' ...
            'state MAT-file'],[1 40],{'SetManually'}));

end
end

% Local function for checking input:
function [wrongInput,input] = checkInput(input)

input = str2double(input{1}); % Converts input to a double.

if isnan(input) || length(input) > 1
    err = errordlg(['The input must be a numeric of length ' ...
        '1 (i.e., a scalar)'],'Error Dialog','modal');
    uiwait(err)
    wrongInput = true;
else
    wrongInput = false; % Correct input was given.
end
end
