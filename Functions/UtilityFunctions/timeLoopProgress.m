% Initialises a progress dialog window for the time loop. 
function timeLoopProgress(start,t,T,status)
%PROGRESSDIALOGWINDOW This function initialises a progress dialog window,
%showing the percentage of completion of the temporal integration. 
%
%   Input: 
%
%   start - A string stating at what time the current simulation was
%   started.
%
%   t - The value of the time variable t at the beginning of the
%   simulation (i.e., the initial time t0). 
%
%   T - The value of the final time in the temporal integration interval. 
%
%   status - Determines if the progress dialog window is enabled or
%   disabled, given as either 'on' or 'off'. 
%
%   Output 
%
%   w - The handle to the progress dialog window; assigned in the callers
%   workspace if status = 'on'. 
%
%   uifig - The handle to the user interface figure window in which the
%   progress dialog window is placed; assigned in the callers workspace if
%   status = 'on'. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

switch status
    case 'on'
        % Creating the user interface figure window in which the progress
        % dialog will be placed:
        uifig = uifigure("WindowStyle","normal","WindowState","normal");

        % Creating the user interface progress dialog with a cancel button:
        w = uiprogressdlg(uifig,'Title','Iterating through time loop', ...
            'Message',[start 'Please wait...'], ...
            'Cancelable','on','ShowPercentage','on');

        % Resizing the the user interface figure window so that it matches
        % the size of the progress dialog:
        uifig.Position(3:4) = [400 185];

        % Adding the value for percentage of completion of the time loop:
        w.Value = t/T;

        % Forcing MATLAB to prioritize finish drawing the progress dialog
        % window:
        drawnow
        pause(.5) % Increase number if the progress dialog doesn't appear.

        % Assigning the outputs to the callers workspace: 
        assignin("caller","w",w)
        assignin("caller","uifig",uifig)
    case 'off'
        % If the progress dialog window is disabled, then the function ends
        % without assigning any outputs to the callers workspace:
        return 
end
end