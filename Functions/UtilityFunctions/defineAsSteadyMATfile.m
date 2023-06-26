function defineAsSteadyMATfile(fileName,Scenario)
%DEFINEASSTEADYMATFILE This function saves a copy of the MAT-file with the
%name fileName in the folder FinalResults (a subfolder of
%SimulationResults), as a MAT-file with the name given by the variable
%version appended with _steady.mat.
%
%   Input:
%
%   fileName - The name of the MAT-file to save as a steady state MAT-file.
%
%   Scenario - The prefix to the name of the steady state MAT-file.
%
%   Output:
%
%   This function has no output. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Creates folders for storing copy steady results, if they do not exist:
mainFolder = 'SimulationResults'; % Name of the main folder.
subFolder = 'FinalResults'; % Name of the subfolder.
if ~isfolder([mainFolder  filesep  subFolder filesep])
    % If there is a folder called SimulationResults, but no subfolder,
    % named FinalResults, the sub folder is
    % created. Furthermore, if neither the main folder nor the
    % subfolder exist, then both are created:
    mkdir([pwd mainFolder  filesep  subFolder filesep],subFolder)
end

destination = [mainFolder filesep subFolder filesep ... % Path to folder
    char(Scenario) '_steady.mat']; % Initial name of MAT-file.

answer = 'No'; % Initiates the while loop.
i = 0; % Initial value for the counter.

% As long as there already exists a file by the name given as input, a
% new name will be made before saving:
while isfile(destination) && string(answer) == "No"
    question = sprintf(['A MAT-file by the name %s already exist, do ' ...
        'you wish to overwrite this file?'],destination);
    answer = questdlg(question,'File already exist','Yes','No','No');
    i = i+1; % Increases the counter by 1. 
    switch answer
        case 'Yes'
            break
        case 'No'
            destination = [mainFolder filesep subFolder filesep ...
                '(' num2str(i) ')' char(Scenario) '_steady.mat'];
    end
end

% Copies the MAT-file to the Results folder and saves it by the name
% defined as destination:
copyfile(correctSlashes(fileName),destination);

% Appends the original file name to the MAT-file: 
originalRun = fileName; 
save(destination,"originalRun","-append")

end

