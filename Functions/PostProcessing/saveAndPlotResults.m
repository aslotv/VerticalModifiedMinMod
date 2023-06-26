% This function saves and plots the results.
function [sol,totalRunTime,RunTime,inputForPlots,Table] = ...
    saveAndPlotResults(initialConditionsInput,Parameters, ...
    Observations,NoDVM,z_0,z_max,dz,u0,un,u,Tol,uMean,TotalPmean, ...
    growth,Production,tspan,time,stopConditionMet,j,Version,Steady, ...
    Scenario,saveResults)
%SAVEANDPLOTRESULTS is a function that saves the results of the simulation
%in a MAT-file, makes the figures from the manuscript and the
%Supporting Information, and gives the option to save them.
%
%   Input:
%
%   initialConditionsInput - The input passed to the function
%   "initialConditions.m".
%
%   Parameters - A structure containing the parameters defined in the main
%   script.
%
%   Observations - A structure containing the observation variables used in
%   the simulation.
%
%   NoDVM - A logical, telling this function whether or not fish was
%   subject to vertical diel migration in this simulation.
%
%   z_0 - The minimum depth of the water column.
%
%   z_max - The maximum depth of the water column.
%
%   dz - The length of the depth steps in the discretised water column.
%
%   u0 - The initial conditions at t = 0.
%
%   un - Solutions from the second to last time step.
%
%   u - Solutions from the last time step.
%
%   Tol - The tolerance for the local error in each time step.
%
%   uMean - The water column average of the solutions.
%
%   TotalPmean - The water column average of the total P.
%
%   growth - The growth in B, A, D, H, C, Z, and F in the last time step.
%
%   Production - The production of B, A, D, H, C, Z and F for each time
%   step in the simulation.
%
%   tspan - The time span for this simulation.
%
%   time - The amount of seconds it took to run the simulation.
%
%   stopConditionMet - Information about why the time loop was terminated.
%
%   j - The sum of consecutive time steps for which the solutions have
%   satisfied the "per-time-step" steady condition.
%
%   Version - Prefix for the name of the MAT-file.
%
%   Steady - A logical telling this function wether or not the steady
%   state condition was met.
%
%   Scenario - A text scalar containing the name of the simulation
%   scenario.
%
%   saveResults - A logical determining whether or not the results are
%   saved in a MAT-file.
%
%   Output:
%
%   sol - A structure with the following fields:
%   * zspan - A vector containing all the depths in the discretised water
%   column.
%   * tspan - The time for each time step.
%   * u0 - The initial values from the first simulation.
%   * uMean - The average water column cocentration for each variable at
%   each time step.
%   * TotalPmean - The average water column concentration of total at each
%   time step.
%   * Production - The production for each osmo- and phagrotrophs from each
%   time step.
%   * u - The solutions from the last time step.
%   * un - The solutions from the second to last time step.
%   * growth - The growth in the osmo- and phagotrophs from the last time
%   step.
%   * steadyCount - The sum of concecutive time steps for which the
%   solutions satisfied the steady tolerance as a function of dt and
%   steadyTime.
%
%   totalRunTime - The amount of time, in seconds, it has taken to run the
%   simulations from the initial time t0 of the first simulation to the
%   final time T of the current simulation.
%
%   RunTime - The amount of time, in seconds, it has taken to run the
%   current simulation.
%
%   inputForPlots - A structure produced by the function
%   "postProcessing.m", containing the input for the plots made by the
%   function "plotResults.m".
%
%   Table - The values for the tables presented in the manuscript,
%   corresponding to the current simulation and the current simulation
%   scenario
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

if strcmp(initialConditionsInput{2},"set")

    % Storing the values z in zSpan, implicitely defined in the
    % pre-time-loop set up:
    sol.zspan = (z_0:dz:z_max)';

    % The time span (from t=t0 to t=T):
    sol.tspan = tspan; 

    % The initial conditions:
    sol.u0 = u0;

    % The mean water column concentration of u at each time t in tspan:
    sol.uMean = uMean; 

    % The mean water column total P concentration at each time t in tspan:
    sol.TotalPmean = TotalPmean;

    % The daily production at each time t in tspan:
    sol.Production = Production; 

    %%% Storing the number of seconds it took to go trough the time loop:
    RunTime = time;
    totalRunTime = time;
else
    % Defining a logical vector in order to find the index of the cell
    % array initialConditions at which the input argument
    % FileName is stored (every element of TF equals false if
    % the value is not stored in the cell array):
    TF = cellfun(@(x) isscalar(string(x)) && ...
        string(x)=="FileName",initialConditionsInput);

    if any(TF)
        % If the optional argument FileName was given as input
        % to the function initialConditions.m, then the file name is
        % defined as the name given as input:
        FileName = correctSlashes(initialConditionsInput{find(TF)+1});
    else
        % If the default file name was used in initialConditions.m, then
        % the file name is defined as the default name:
        FileName = 'SimulationResults\FinalResults\Baseline_steady.mat';
    end

    % Loading the contents of the MAT-file and storing it in a structure:
    loadedVariables = load(correctSlashes(FileName));

    % Extracting the solution structure from the loadedVariables structure:
    sol = loadedVariables.sol;

    % Defining a logical vector in order to find the index of the cell
    % array initialConditions at which the input argument
    % ContinuedSimulation is stored (every element of TF equals false if
    % the value is not stored in the cell array):
    TF = cellfun(@(x) isscalar(string(x)) && ...
        string(x)=="ContinuedSimulation",initialConditionsInput);

    % Checks if the current simulation was a new simulation or a continued
    % simulation:
    if any(TF)
        % If the optional argument ContinuedSimulation was given as input
        % to the function initialConditions.m (or added to the cell array
        % initialConditionsInput within the function due to missing fields
        % in sol), then the logical newSimulation is defined as the
        % negation of ContinuedSimulation:
        newSimulation = ~initialConditionsInput{find(TF)+1};
    else
        % If ContinuedSimulation was not given as input to
        % initialConditions.m (and not added within the funciton), then
        % newSimulation is defined as the negation of the default value of
        % ContinuedSimulation (ContinuedSimulation = true):
        newSimulation = false;
    end

    if newSimulation

        % Storing the values z in zSpan, implicitely defined in the
        % pre-time-loop set up:
        sol.zspan = (z_0:dz:z_max)';

        % The time span (from t=t0 to t=T):
        sol.tspan = tspan;

        % The initial conditions:
        sol.u0 = u0;

        % The mean water column concentration of u at each time t in tspan:
        sol.uMean = uMean;

        % The mean water column total P concentration at each time t in
        % tspan:
        sol.TotalPmean = TotalPmean;

        % The daily production at each time t in tspan:
        sol.Production = Production; 

        %%% Storing the number of seconds it took to go trough the time
        %%% loop:
        RunTime = time;
        totalRunTime = time;

    else

        % Extracts the variable totalRunTime from the loadedVariables
        % structure if it exists or defines it as zero if not:
        if isfield(loadedVariables,'totalRunTime')
            totalRunTime = loadedVariables.totalRunTime;
        else
            totalRunTime = 0;
        end

        % The time span (from t=T_old to t=T):
        sol.tspan = [sol.tspan tspan(2:end)+sol.tspan(end)];

        % The mean water column concentration of u at each time in tspan:
        sol.uMean = [sol.uMean; uMean(2:end,:)];

        % The mean water column total P concentration at each time t in
        % tspan:
        sol.TotalPmean = [sol.TotalPmean; TotalPmean(2:end,:)];

        % The production at each time t in tspan:
        sol.Production = [sol.Production; Production(2:end,:)];

        %%% Storing the time it took to iterate through the time loop:
        RunTime = time; % The time this run took.
        totalRunTime = totalRunTime + time; % The combined run time.
    end
end

% Storing the results from final time step T:
sol.u = u;
sol.un = un;
sol.growth = growth;

% Storing the sum of consecutive time steps with solutions meeting the
% steady tolerance:
sol.steadyCount = j; % j is in hours.

% Using the function "postProcessing.m" to define and calculate the input
% data for plotting and for the tables in the manuscript:
[inputForPlots,Table] = postProcessing('Solutions',sol, ...
    'Parameters',Parameters,'Observations',Observations,'NoDVM',NoDVM);

if saveResults
    % The current date and time:
    thisYear = string(datetime('now','format','(yyyy)'));
    todaysDate = string(datetime('now','format','MMMMd'));
    timeNow = string(datetime('now','format','h.mma'));

    % Creates folders for storing simulations, if they do not exist:
    mainFolder = "SimulationResults"; % Name of the main folder.
    subFolder = todaysDate + thisYear; % Name of the subfolder. 
    if ~isfolder(mainFolder + filesep + subFolder)
        % If there is a folder called SimulationResults, but no subfolder,
        % named with today's date and the current year, the sub folder is
        % created. Furthermore, if neither the main folder nor the
        % subfolder exist, then both are created: 
        mkdir(pwd + string(filesep) + mainFolder + string(filesep), ...
            subFolder)
    end

    % The parameters and observations, the pre-time-loop set up, and the
    % solutions are saved in the subfolder named with today's date and the
    % current year. The mat-file is named with the current time (e.g.,
    % 1.01PM):
    filePath = mainFolder + string(filesep) + subFolder + string(filesep);
    fileName = Version + timeNow;

    % Saving the results:
    save(filePath + fileName + ".mat", 'stopConditionMet', ...
        'initialConditionsInput','Parameters', 'Observations','sol', ...
        'totalRunTime','RunTime','inputForPlots','Table','Scenario','Tol')

    % Print message with information on where the results have been saved:
    fprintf("The results of this run have been saved in the file %s" ...
        + ", in the sub folder %s, in the folder ""%s"".\n", ...
        fileName + ".mat",subFolder,mainFolder) %#ok<PRTCAL>
end

% Gives the option to save a cpopy of the results in a MAT-file in the
% folder "FinalResults" with the prefix given as the input Scenario and the
% postfix "_steady.mat":
if Steady
    answer = questdlg(['Save a copy of the results in a MAT-file in the'...
        ' folder FinalResults?'],'Save steady state results', ...
        'Yes','No','Yes');
    if answer == "Yes"
        defineAsSteadyMATfile(filePath + fileName + ".mat",Scenario)
    end
end

% Producing the figures in the manuscript and in the Supporting
% Information:
figures = plotResults('InputForPlots',inputForPlots);

% Gives an option to save the figures:
answer = questdlg('Do you want to save the figures?', ...
    'Option to save figures','Yes','No','Yes');

switch answer
    case 'Yes'
        mkdir(filePath,"Figures_" + timeNow) % Creates folder.
        for i = 1:length(fieldnames(figures))
            savefig(filePath + "Figures_" + timeNow + ...
                string(filesep) + string(get(gcf,'Name')) + ".fig")
            close(gcf)
        end
        fprintf(['\nThe figures are saved in the same place as the ' ...
            'results, in a folder with the name %s.\n'], ...
            "Figures_" + timeNow)

        % Gives the option to save a copy of the figures in the folder
        % "FinalResults" with the prefix "Figures_", the infix given as the
        % input Scenario and the postfix "_steady":
        if Steady
            answer = questdlg(['Save a copy of the figures in the'...
                ' folder FinalResults?'],['Save copy of steady state ' ...
                'figures'],'Yes','No','Yes');
            if answer == "Yes"
                defineAsSteadyFigures(filePath + "Figures_" + timeNow, ...
                    Scenario)
            end
        end

    case 'No'
        fprintf(['\nFigures were not saved; use function ' ...
            'plotResults.m to make new figures and save them.\n'])
        help plotResults
end
end

