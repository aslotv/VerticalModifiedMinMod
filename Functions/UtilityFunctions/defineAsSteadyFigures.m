function defineAsSteadyFigures(folderName,Scenario)
%DEFINEASSTEADYFIGURES This function saves a copy of the figures in the
%folder whos name is given by the input argument folderName, in a folder
%named by the prefix "Figures_", the infix given by the input argument
%Scenario, and the postfix "_steady", in the folder FinalResults (a
%subfolder of SimulationResults).
%
%   Input: 
%
%   folderName - The name of the folder in which the figures to be copied
%   are stored.
%
%   Scenario - The name of the simulation scenario for which the figures
%   were produced. 
%
%   Output: 
%
%   This function has no output. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 


% Creates folders for storing copy of steady results, if they do not exist:
mainFolder = 'SimulationResults'; % Name of the main folder.
subFolder = 'FinalResults'; % Name of the subfolder.
if ~isfolder([mainFolder  filesep  subFolder filesep])
    % If there is a folder called SimulationResults, but no subfolder,
    % named FinalResults, the sub folder is
    % created. Furthermore, if neither the main folder nor the
    % subfolder exist, then both are created:
    mkdir([pwd filesep mainFolder filesep],subFolder)
end

destination = [mainFolder filesep subFolder filesep ... % Path to folder
    'Figures_' char(Scenario) '_steady']; 

answer = 'No'; % Initiates the while loop.
i = 0; % Initial value for the counter.

% As long as there already exists a folder by the name given as input, a
% new name will be made before saving:
while isfolder(destination) && string(answer) == "No"
    question = sprintf(['A folder by the name %s already exist, do ' ...
        'you wish to overwrite this folder?'],destination);
    answer = questdlg(question,'Folder already exist','Yes','No','No');
    i = i+1; % Increases the counter by 1. 
    switch answer
        case 'Yes'
            break
        case 'No'
            destination = [mainFolder filesep subFolder filesep ...
                '(' num2str(i) ')' 'Figures_' char(Scenario) '_steady'];
    end
end

% Copies the figure folder to the Results folder and saves it by the name
% defined as destination:
copyfile(correctSlashes(folderName),destination);

end

