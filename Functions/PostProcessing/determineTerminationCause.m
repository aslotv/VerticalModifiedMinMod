function [stopConditionMet,saveResults] ...
    = determineTerminationCause(j,steadyTime,steadyTol,tSteady,T,t,dt, ...
    abortIntegration,dtThresh,nonNegTol,e)
%DETERMINETERMINATIONCAUSE This function determines why the temporal
%integration was terminated and stores it, as well as giving the option to
%save the results of the simulation in a MAT-file. Furthermore, this
%function plays a sound, alerting to the fact that the temporal integration
%has been terminated; the type of sound played, depends on the reason for
%termination.
%
%   Input:
%
%   j - The sum of consecutive time steps for which the solutions have
%   satisfied the "per-time-step" steady condition.
%
%   steadyTime - The desired length of time for which the solutions'
%   temporal variation must be within a desired tolerance interval.
%
%   steadyTol - A "per-time-step" steady tolerance, defined as a function
%   of steadyTime and the current time step.
%
%   tSteady - The time t at which the solutions reached steady state.
%
%   T - The final time of the temporal integration interval.
%
%   t - The value of the time t when the temporal integration was
%   terminated. 
%
%   dt - The length of the last successful time step. 
%
%   abortIntegration - A logical telling this function whether or not the
%   temporal integration was aborted due to the length of the time steps
%   falling below the threshold. 
%
%   dtThresh - The threshold for the length of the time steps as a function
%   of the current time t. 
%
%   nonNegTol - An error tolerance, which equals Tol (the tolerance for the
%   local error) if all non-negative variables have only non-negative
%   values, and equals zero if any of the non-negative variables have any
%   negative values. 
%
%   e - The error in the current time step, as a function of the local
%   error and the "negative" error. 
%
%   Output: 
%
%   stopConditionMet - A text scalar containing information about why the
%   temporal integration was terminated. 
%
%   saveResults - A logical telling the function saveAndPlotResults.m
%   whether or not to save the results from the simulation as a MAT-file. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

if j >= steadyTime

    % If steady state was reached:
    stopConditionMet = sprintf("Reached steady state with " + ...
        "steadyTol = %s, at t = %s. Stopped the time loop at " + ...
        "T = hours(days(%s)).", ...
        string(steadyTol(dt)),string(tSteady),string(days(hours(T))));
    stopSound = 'train';

elseif t == T

    % If final T was reached before steady state:
    stopConditionMet = sprintf("Reached final time " + ...
        "T = hours(days(%s))",string(days(hours(T))));
    stopSound = 'train';

elseif abortIntegration

    % If integration was aborted:
    stopConditionMet = "Simulation was aborted; unable to " + ...
        "meet tolerance without dt becoming too small.";
    stopSound = 'gong';

    % Issuing a warning if the integration was aborted:
    warning(['The integration is aborted at t = %g, due to ' ...
        'the value of the time step falling below the ' ...
        'threshold (dtThres = %g) when trying to meet the ' ...
        'error tolerance Tol = %g (the error from the last ' ...
        'attempt was %g).'],t,dtThresh(t),nonNegTol,e(n))

    % Gives the option to save the results from the aborted integration:
    saveResults = questdlg(['Do you want to save the results ' ...
        'from the aborted simulation (from t0 to the last whole day' ...
        ' for which a solution was successfully found)?'], ...
        'Save results','Yes','No','Yes');

else
    % If simulation was canceled:
    stopConditionMet = "Simulation was canceled.";
    stopSound = 'gong';

    % Gives the option to save the results from the canceled integration:
    saveResults = questdlg(['Do you want to save the results ' ...
        'from the canceled simulation?'], ...
        'Save results','Yes','No','Yes');
end

% Prints information about why the time loop ended:
fprintf("\n" + stopConditionMet + "\n")

try % Plays a sound when the time loop has ended.
    playSound(stopSound)
    pause(10)
catch ME
    warning(ME.identifier,'%s',ME.message)
end

% If the simulation was not canceled or aborted the option to save the
% results is also given. However, if neither option is chosen before the
% timer runs out, the results are saved automatically by the function
% "saveAndPlotResults.m":
if ~exist('saveResults','var')

    % Defining the number of seconds for which it is still possible to opt
    % out of saving:
    timeOut = 30; 

    % Creating the user interface figure window in which the confirmation
    % dialog window will be placed:
    uifig = uifigure("WindowStyle","normal","WindowState","normal");

    % Resizing the the user interface figure window so that it matches the
    % size of the confirmation dialog window:
    uifig.Position(3:4) = [400 150];

    % Creating the user interface progress dialog with a cancel button:
    msg = @(timeOut,timeElapsed) sprintf(['Saving the results in %g ' ...
        'seconds, press cancel to abort saving.'], ...
        floor(timeOut-timeElapsed)); 

    w = uiprogressdlg(uifig,'Title','Do you want to save the results?', ...
        'Message',msg(timeOut,0), ...
        'Cancelable','on');

    w.Value = 0; % Sets initial value for the progress bar. 

    % Forces MATLAB(R) to finish drawing progress dialog before continuing:
    drawnow
    pause(.5)

    startedTimer = tic; % Starts the timer. 

    while timeOut-toc(startedTimer)>0
        % Updates the uiprogressdlg:
        w.Message = msg(timeOut,toc(startedTimer)); 
        w.Value = min(1,toc(startedTimer)/timeOut);

        % Stops while loop if run is canceled:
        if w.CancelRequested 
            saveResults = false; 
            break
        end
    end
    stopped = toc; %#ok<NASGU> 
    delete(uifig)

    % If the cancel button was not pressed, then the results will be saved:
    if ~exist('saveResults','var'); saveResults = true; end
end
end