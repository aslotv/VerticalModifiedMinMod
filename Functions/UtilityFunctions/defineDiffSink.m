% Defines the diffusivity and sinking vectors. 
function [Kappa,Nu] = defineDiffSink(kappa,nu,N_vars)
%DEFINEDIFFSINK is a function that puts kappa and nu in vector form, in
%order to simplify the calculation of the right hand side of the partial
%differential equations in the Vertical Modified MinMod. 
%
%   Input: 
%
%   kappa - The turbulent diffusivity. 
%
%   nu - A structure containing the non-zero sinkning velocities. 
%
%   N_vars - The number of state variables or equations in the system of
%   partial differential equations of the Vertical Modified MinMod. 
%
%   Output: 
%
%   Kappa - The diffusivity on vector form. 
%
%   Nu - The sinking on vector form. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% The diffusivity: 
assert(isscalar(kappa),['The function defineDiffSink.m is' ...
                        ' does not accept kappa varying with depth,' ...
                        ' nor different values of kappa for different' ...
                        ' state variables yet.'])
Kappa = ones(1,N_vars)*kappa;  % kappa assumed constant. 
Kappa(7) = 0; % No diffusion in fish equation. 

% The sinking: Only slow and fast sinking detritus, and particulate organic
% silicate has nonzero sinking coefficient:
Nu = zeros(1,N_vars); % Storage vector (for increased index control).
Nu(8:9) = [nu.Det_s nu.Det_f]; % For slow and fast sinking detritus.
Nu(14) = nu.Det_f; % S_opal has the same sinking rate as Det_f. 
end