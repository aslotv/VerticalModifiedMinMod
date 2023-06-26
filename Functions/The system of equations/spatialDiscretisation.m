function [N_z,dz,zspan,zrange] = spatialDiscretisation(z_0,z_max,dz)
%SPATIALDISCRETISATION is a function that defines the spatial
%discretisation of the system of PDEs based on the minimum and maximum
%depth of the water column, and the desired length of the depth steps. 
%   
%   Input: 
%
%   z_0 - The minimum depth of the water column.
%
%   z_max - The maximum depth of the water column.
%
%   dz - The desired length of the depth steps. 
%
%   Output: 
%
%   N_z - The number of depths in the discretised water column. 
%
%   dz - The length of the depth steps in the water column. Note that the
%   value of dz is changed within this function if the value given as input
%   results in failure to define a depth span including the maximum depth
%   of the water column. A warning will be issued if this is the case. 
%
%   zspan - A vector of all the depths in the discretised water column. 
%
%   zrange - The depth range of the water column. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

zrange = z_max - z_0; % The depth range in the water column. 

N_z = zrange/dz + 1; % Number of depths in zspan.

if ceil(N_z) ~= N_z % Redefines dz and N_z if N_z is not an integer.

    dz = zrange/ceil(N_z); % Defining dz so that N_z is an integer.

    warning('N_z=%.2f is not an integer, using altered value of dz ~ %g', ...
        N_z,dz);

    N_z = ceil(N_z) + 1;
end

zspan = (z_0:dz:z_max)'; % A vector with the discretised depths.

end

