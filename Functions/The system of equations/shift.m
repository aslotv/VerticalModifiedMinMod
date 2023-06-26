% Copies shifted elements in un while handeling boundary conditions.
function [u_shiftUpDiffusion,u_shiftDownDiffusion,u_shiftDownSinking] ...
             = shift(un,s,N_vars,P,S_b,kappa)
%SHIFT is a function that shifts the elements in un one step up and one
%step down, so that u_shiftUp(i,:)=un(i+1,:) and
%u_shiftDown(i,:)=un(i-1,:). In addition, the boundary conditions are
%accounted for in the function defineBoundaryConditions.m, and added to the
%shifted matrices.
%
%   Input: 
%
%   un - The solutions from the previous time step. 
%
%   s - The net growth in the state variables from the current time step. 
%
%   N_vars - The number of state variables or equations in the system of
%   partial differential equations of the Vertical Modified MinMod. 
%
%   P - A structure containing the photosynthetic quotient, inorganic
%   phosphate concentration below deep boundary (z_max), and the
%   particulate losses.
%
%   S_b - The inorganic silicate concentration below z = z_max. 
%
%   kappa - The turbulent diffusivity for all the state variables. 
%
%   Output: 
%
%   u_shiftUpDiffusion - The first term in the central difference
%   approximating the diffusion term in the system of PDEs, altered to
%   satisfy boundary conditions.
%
%   u_shiftDownDiffusion - The last term in the central difference
%   approximating the diffusion term in the system of PDEs, altered to
%   satisfy boundary conditions.
%
%   u_shiftDownSinking - The last term in the backward difference
%   approximating the sinking term in the system of PDEs, altered to
%   satisfy boundary conditions.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Using the function defineBoundaryConditions to account for boundary
% conditions in the shifted matrices:
[diffusionSurface,diffusionBottom,sinkingSurface,sinkingBottom] = ...
    defineBoundaryConditions(un,s,kappa,P,S_b,N_vars);

% The first term in the central difference approximating the diffusion term
% in the system of PDEs, is defined as u_shiftUpDiffusion. The last row of
% this matrix, is defined in such a way that the diffusion at the bottom
% boundary, is in line with the boundary conditions (see
% defineBoundaryConditions.m for details). For the remaining rows,
% u_shiftUpDiffusion(i,:) = un(i+1,:), with i=1,2,...,N_z-1.
u_shiftUpDiffusion = [un(2:end,:); diffusionBottom];

% The first term in the central difference approximating the diffusion term
% in the system of PDEs, is defined as u_shiftUpDiffusion. The first row of
% this matrix, is defined in such a way that the diffusion at the surface
% boundary, is in line with the boundary conditions (see
% defineBoundaryConditions.m for details). For the remaining rows,
% u_shiftDownDiffusion(i,:) = un(i-1,:), with i=2,...,N_z.
u_shiftDownDiffusion = [diffusionSurface; un(1:end-1,:)];

% The last term in the backward difference approximating the sinking term
% in the system of PDEs, is defined as u_shiftDownSinking. The first and
% last row of this matrix, are defined in such a way that the sinking
% through the boundary, is in line with the boundary conditions (see
% defineBoundaryConditions.m for details). For the remaining rows,
% u_shiftDownSinking(i,:) = un(i-1,:), with 1<i<N_z.
u_shiftDownSinking= [sinkingSurface; un(1:end-2,:); sinkingBottom];

end
