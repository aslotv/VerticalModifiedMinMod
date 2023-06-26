%% Defining the Observation MAT-files Used in the Simulations
% There are three MAT-files used in the simulations. These three files
% contains observed diel vertical fish distribution, and vertical profiles
% of temperature and oxygen concentrations. The variables in these files
% need some minor alterations before they can be used in the simulations.
% This script contains the code used for said alterations, as well as
% saving the altered variables in new MAT-files that can be used directly
% in the script "VerticalModifiedMinMod_dla_asl.m".
%
%   Authors: Anita Stene Loetvedt and Dag L. Aksnes.

%% The observed temperature:
% The altered temperature profile is only used in defining parameters.

% The temperature profile from 1m to 700m: 
load Observations\FORCE_temperature_700m.mat temp

% Setting temperature at z=0 equal to temp(1):
if isrow(temp); temp = temp'; end
temp = [temp(1); temp];

% Saving the altered temperature profile: 
save('Observations\temperatureProfile_0to700m.mat','temp')

%% The observed oxygen concentrations:
% The altered oxygen profile is only used in defining initial conditions. 

% The oxygen profile from 1m to 700m: 
load Observations\OBS_den_chl_oxy_pen.mat oxy 

% At z = 0, initial conditions for oxygen is set to 100%
% saturation at temperature 27.7 and salinity 39.3:
if isrow(oxy); oxy = oxy'; end 
oxy = [4.4; oxy]*(44.7e3); % From ml/L to nM. 

% Saving the altered oxygen profile: 
save('Observations\oxygenProfile_0to700m.mat','oxy'); 

%% The observed diel vertical fish distribution:
% The altered diel vertical fish distribution is used in defining initial
% conditions, in defining the right hand side of the system of semi
% discretised PDEs, and in the post processing. 

% The average 24h relative vertical fish distribution:
load Observations\fishdistribution_682m fishdist % Load distribution.

% All missing elements in the vector containing the vertical distribution
% are equal to zero:
if isrow(fishdist); fishdist = fishdist'; end
fishDistribution = [0; fishdist; zeros(18,1)];
clear fishdist % Cleans up workspace.

% Saving the fish distribution from 0m to 700m:
save('Observations\fishDistribution_0to700m.mat',"fishDistribution");
