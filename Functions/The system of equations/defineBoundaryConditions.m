function [diffusionSurface,diffusionBottom, ...
    sinkingSurface,sinkingBottom] = ...
    defineBoundaryConditions(un,s,kappa,P,S_b,N_vars)
%DEFINEBOUNDARYCONDITIONS This function defines the boundary conditions for
%the system of partial equations in the Vertical Modified MinMod, by
%defining row vectors that will replace the top and bottom rows of the
%shifted matrices used to estimate the diffusion and convection terms in
%the system of PDEs. 
%
%   Input: 
%
%   un - The solution from the previous time step. 
%
%   s - The net growth from the current time step. 
%
%   kappa - The turbulent diffusivity. 
%
%   P - A structure containing the photosynthetic quotient, inorganic
%   phosphate concentration below deep boundary (z_max), and the
%   particulate losses.
%
%   S_b - The inorganic silicate concentration below z = z_max. 
%
%   N_vars - The number of equations in the system of PDEs. 
%
%   Output : 
%   
%   diffusionSurface - The values of u(t_n) at z = z_0 - dz, defined in
%   such a way that boundary conditions regarding diffusion at upper
%   boundary are satisfied.  
%
%   diffusionBottom - The values of u(t_n) at z = z_max + dz, defined in
%   such a way that boundary conditions regarding diffusion at bottom
%   boundary are satisfied. 
%
%   sinkingSurface - The values of u(t_n) at z = z_0 - dz, defined in
%   such a way that boundary conditions regarding sinking through upper
%   boundary are satisfied.
%
%   sinkingBottom- The values of u(t_n) at z = z_max - dz, defined in
%   such a way that boundary conditions regarding sinking through bottom
%   boundary are satisfied.
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Setting diffusion from above surface boundary equal to diffusion from
% surface boundary up through said boundary:
diffusionSurface = un(1,:); 

% For oxygen (O), the diffusion at the upper boundary is defined so that
% the oxygen concentration at zero depth is constant with respect to time
% (Dirichlet boundary condition u15(t,z_0) = u15(t0,z_0) = 100% oxygen
% saturation):
diffusionSurface(15) = 2*un(1,15) - un(2,15) - s(1,15)/kappa(15);

% Setting diffusion from below bottom boundary equal to diffusion from
% bottom boundary down through said boundary:
diffusionBottom = un(end,:); 

% For dissolved inorganic phosphate (P) and silicate (S), there is
% a constant supply from below z_max:
diffusionBottom(:,11) = 1000*P.b; % Sets concentration for P below z_max.
diffusionBottom(:,13) = 1000*S_b; % Sets concentration for S below z_max.

% Setting sinking from above surface boundary equal to zero: 
sinkingSurface = zeros(1,N_vars);

% No constraints on sinking through bottom boundary: 
sinkingBottom = un(end-1,:); 
end
