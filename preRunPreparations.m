function FolderName = preRunPreparations(WorkspaceVars,FolderName,options)
%PRERUNPREPARATIONS is a function that gives the option to save the
%variables in the callers workspace, if any, and adds the folders
%(including the subfolders) given as input to the MATLAB(R) search path.
%
%   Input: 
%
%   WorkspaceVars - One of two MATLAB(R) functions, who or whos; if either
%   of these are defined as a variable before passing it to this function,
%   said variable will be saved in the MAT-file along with any pre-existing
%   variables.
%
%   FolderName - This is a so called repeating argument, which means that
%   names of folders to add to the MATLAB(R) search path can be given as a
%   comma separated list. There are two things to note her: The first is
%   that if the folder is not in the current directory, the folder name
%   must include the file path. The second is that this function will also
%   add any sub folders in the folders (with the folder name) given as
%   input to the MATLAB(R) search path, not just the parent folder. 
%
%   Optional input: 
%
%   RunAsExample (RunAsExample = false, as default) - Set the value of
%   RunAsExample equal to true if you wish to run this function without
%   your workspace variables being saved and without making any alterations
%   to the MATLAB(R) search path.
%
%   Output: 
%
%   FolderName - A cell containing the names of the folders added to the
%   MATLAB(R) search path. 
%
%   Authors: Anita Stene Loetvedt and Dag L. Aksnes.

arguments
    WorkspaceVars {mustBeA(WorkspaceVars,{'cell','struct'})}
end

arguments(Repeating)
    FolderName {mustBeFolder}
end

arguments
    options.RunAsExample (1,1) {mustBeNumericOrLogical} = false; 
end

%% Giving option to save variables in callers workspace (if any)
if ~isempty(WorkspaceVars) % If there are any variables in the workspace.
    % Gives an option to save the workspace:
    answer = questdlg('Do you want to save your workspace?', ...
        'Option to save workspace','Yes','No','Yes');
    switch answer
        case 'Yes'
            if options.RunAsExample
                % If the function is run as an example, display a message
                % explaining what would have happened when opting to save: 
                msg = sprintf(['Your workspace variables would have ' ...
                    'been saved in a MAT-file named "%s".'], ...
                    "YourWorkspace" + ...
                    string(datetime('now','format','dMMMh.mma')) ...
                    + ".mat");
                uiwait(msgbox(msg,'Example','help','modal'));
            else
                % If the function is not run as an example, save the base
                % workspace variables in a MAT-file, named "YourWorkspace"
                % with the current date and time as a postfix: 
                evalin("base", ...
                    "save(string('YourWorkspace') + " + ...
                    "string(datetime('now','format','h.mm.a')) + " + ...
                    "string('.mat'))")
            end
        case 'No'
            if options.RunAsExample
                fprintf(['\nThe following variables would have been ' ...
                    'cleared from memory without saving them in a ' ...
                    'MAT-file first:\n'])
                evalin("base","whos")
            else
                fprintf('\nWorkspace was not saved.\n')
            end
    end
end

%% Adding folders from current directory to MATLAB(R) search path:

if options.RunAsExample % If the function is run as an example. 
    return % Ends function without altering the MATLAB(R) search path. 
end

for i = 1:length(FolderName)
    
    % Creating a full path to the folder, including possible sub folders:
    FolderPath = genpath([pwd filesep char(FolderName{i})]);

    % Adding the folder whose full path is saved in FolderPath, including
    % possible sub folders, to the MATLAB(R) search path, if it's not
    % already on it:
    if ~contains(path,caseSensitivePattern(FolderPath))
        addpath(FolderPath)
    end
end
end