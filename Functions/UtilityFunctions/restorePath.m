function restorePath(FolderName)
%RESTOREPATH This function gives the option to remove the folders whose
%name is given as input, from the MATLAB(R) search path. 
%
%   Input: 
%
%   FolderName - A cell array where each element is the name of a folder
%   previously added to the MATLAB(R) search path. For each of the folders
%   whose name is in this cell array, the option to remove the folder from
%   the search path is given.
%
%   Output: 
%
%   This function has no output. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Creating the user interface figure window in which the confirmation
% dialog window will be placed:
uifig = uifigure("WindowStyle","normal","WindowState","normal");

% Resizing the the user interface figure window so that it matches the size
% of the confirmation dialog window:
uifig.Position(3:4) = [400 160];

for i = 1:length(FolderName)
    % Creating the user interface confirmation dialog window:
    msg = ['Remove the folder ' char(FolderName{i}) ' and its sub ' ...
        'folders from MATLAB(R) search path?'];
    answer = uiconfirm(uifig,msg,['Comfirm remove folder from ' ...
        'search path'],'Options',{'Remove','Keep'}, ...
        'DefaultOption','Remove','CancelOption','Keep');

    switch answer
        case 'Remove'
            % Removing folders added to the MATLAB search path:
            rmpath(genpath([pwd filesep char(FolderName{i})]))
            uialert(uifig,['The folder named "' char(FolderName{i}) '"' ...
                ' was removed from the MATLAB search path '], ...
                'MATLAB search path','Modal',true,'Icon','info', ...
                'CloseFcn',@resume);
            uiwait(uifig)
        case 'Keep'
            uialert(uifig,['The folder named "' ...
                char(FolderName{i}) '" remains in the MATLAB search ' ...
                'path until removed manually or until you quit MATLAB'],...
                'MATLAB search path','Modal',true,'Icon','info', ...
                'CloseFcn',@resume);
            uiwait(uifig)
    end
    clear answer msg
end
delete(uifig)

    function resume(src,event)
        uiresume(uifig)
    end
end