function [filePath] = correctSlashes(filePath)
%CORRECTSLASHES is a function that replaces backslashes (\) with forward
%slashes (/), in the given file path (filePath), if running function on Mac
%OS X or UNIX version of MATLAB, and replaces forward slashes (/) with
%backslashes (\) in the given file path if running function on PC (Windows)
%version of MATLAB. 
%
%Note: Using this function on a filePath that needs no corrections will
%only result in the function returning a file path that is equivalent to
%the one given as input. Hence, this function can be used to ensure correct
%file separators in file paths. 
%
%   Input: filePath - The file path for which the slashes needs (or may
%   need) to be corrected.
%
%   Output: 
%   filePath - The file path with corrected slashes. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

slashes = {'\','/'}; 

filePath = replace(filePath,slashes,filesep);

end

