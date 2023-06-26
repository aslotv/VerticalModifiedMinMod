% Calculates the values in the tables and the input for the plots.
function [inputForPlots,Table] = postProcessing(optional)
%POSTPROCESSING is a function that calculates the values presented in the
%tables in the manuscript, as well as defining and calculating data for
%plotting.
%
%   Inputs:  
%   All inputs are optional, however, if any of the optional inputs
%   Solutions, Parameters, and Observations are given, then all three must
%   be given. Give these optional inputs if you don't want to load them
%   from a MAT-file. 
%
%   FileName (FileName  = 'Baseline_steady.mat' as default) - The name of
%   the MAT-file where the results from the simulations are saved.
%
%   Solutions - A structure containing the following fields: 
%   * tspan - The times t in tspan.
%   * zspan - All the depths in the water column. 
%   * u - The solutions from the last time step.
%   * un - The solutions from the second to last time step.
%   * growth  -  The growth in the osmo- and phagotrophs from the last time
%   step.
%   * Production - The production for the osmo- and phagotrofs for each
%   t in tspan.
%   * uMean - The mean concentrations from each time t in tspan.
%   * TotalPmean - The mean concentrations of total P for each time t in
%   tspan.
%
%   Parameters - A structure containing the following fields:
%   * alfa - A structure containing the affinities and clearance rates.
%   * sigma - Mesozooplankton selectivity factor for C relative to D.
%   * mu - A structure containing the maximum growth rates for B, A, D, H,
%   C, Z, and F, as well as storage vectors for the growth rates for B, A,
%   and D.
%   * Y - A structure containing the yields of H, C, Z, F, and B.
%   * I- A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as storage vectors for the ingestion rates and the
%   denominator of the fractions defining them.
%   * H - A structure containing the half saturation constants.
%   * limFrac - A structure containing storage vectors for the limitaion
%   fractions.
%   * limDepths - A structure containing logical storage vectors for the
%   logical depth indices corresponding to the depths at which each of the
%   limiting factors apply.
%   * f - A structure containing photosynthetic carbon overflow factor, and
%   the fraction of total loss that goes to detritus and that is respired.
%   * delta - A structure containing the mortality rates for D and F, and
%   for Z if NoFish = true.
%   * k - A structure containing the dissolution, leakage and fragmentation
%   rates, as well as empirical coefficients related to chlorophyll light
%   absorption.
%   * nu - A vector containing the sinking velocities for all variables
%   (both zero and non-zero).
%   * K - A structure containing the light attenuation due to clear water
%   and do to other factors than clear water and chlorophyll concentration.
%   * P - A structure containing the photosynthetic quotient, inorganic
%   phosphate concentration at deep boundary (z_max), and the particulate
%   losses.
%   * R - A structure containing the respiratory quotient, and respiration
%   losses.
%   * rho - A structure containing the stoichiometric ratios and conversion
%   factors.
%   * E_0 - The surface irradiance. 
%   * kappa - A vector containing the turbulent diffusivity for all
%   variables (both zero and non-zero).
%   * S_b - The inorganic silicate concentration at z = z_max.
%
%   Observations - A structure containing the following fields:
%   * fishDistribution - The relative average diel distribution of fish in
%   the water column.
%   * oxygenProfile - The vertical oxygen profile in the water column. 
%   * temperatureProfile - The vertical temperature profile in the water
%   column.
%
%   noDVM (NoDVM = false as default) - Tells the function whether or not
%   the simulation was done with fish subjected to diel vertical migration.
%
%   Output:
%
%   inputForPlots - A structure containing the input for the plots made by
%   the function "plotResults.m", which has the following fields:
%   
%   * meanConcentrations - The input for the average water column
%   concentrations figure. 
%
%   * verticalProfiles - The input for the vertical profiles of the
%   solutions figure. 
%
%   * carbonFluxAndFishProcesses - The input for the carbon flux and fish
%   processes figure. 
%
%   * relativeCarbonFluxes - The input for the relative carbon fluxes
%   figure. 
%
%   * ObservationsAndSimulations - The input for the observations and
%   simulations figure. 
%   
%   Table - A table containing the values from the current run which are
%   presented in the tables in the manuscript. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%% Arguments block
arguments
    optional.FileName (1,:) {mustBeText} = 'Baseline_steady.mat'
    optional.Solutions {mustBeA(optional.Solutions,'struct')}
    optional.Parameters {mustBeA(optional.Parameters,'struct')}
    optional.Observations {mustBeA(optional.Observations,'struct')}
    optional.NoDVM {mustBeNumericOrLogical} = false;
end

%%% Asserting that necessary input is given if results are passed as input

% The complete set of inputs which must be given as input, if any of its
% members are given:
InputSet = {'Solutions','Parameters','Observations'};

% Checking if the complete set of inputs are given, if one of its
% members are given:
AnyInputInSetGiven = any(cellfun(@(x) isfield(optional,x),InputSet));
isCompleteInputsetGiven = all(cellfun(@(x) isfield(optional,x),InputSet));

% Throws en error if one or more of the inputs in InputSet are given, but
% not all:
if AnyInputInSetGiven && ~isCompleteInputsetGiven && ~optional.NoDVM

    errID = 'MATLAB:narginchk:notEnoughInputs';
    errMsg = sprintf(['If any of the arguments ''Solutions'', ' ...
        '''Parameters'', or ''Observations'' are given as input, ' ...
        'then all of them must be given as inputs:\n' ...
        '[inputForPlots,Table] = postProcessing(''Solutions'',sol,' ...
        '''Parameters'',Parameters,''Observations'',Observations)']);

    ME = MException(errID,errMsg);
    throw(ME)
end

%% Extracting variables from MAT-file/optional input structure

% Solutions from t=0, to t = T, the parameters, and the observations:
if ~AnyInputInSetGiven
    % The simulation results saved in the MAT-file with the name FileName,
    % are loaded into the workspace for post processing, if they are not
    % given as input:
    load(optional.FileName,'sol','Parameters','Observations')
else
    sol = optional.Solutions;
    Parameters = optional.Parameters;
    Observations = optional.Observations;
end

%% Defining size solutions matrix and length of depth step
% Number of depths in zspan, N_z, and number of variables, N_vars:
[N_z,N_vars] = size(sol.u);

% Length of depth step: 
dz = diff(sol.zspan); 

% If all depth steps are equal, define dz as a scalar:
if min(dz) == max(dz); dz = dz(1); end

%% Production and fish biomass at the last time step

% Primary daily production of autotrophs (in mg C m^(-2) d^(-1)):
NPP = 24*sum(sol.Production(end,2:3),2).*Parameters.rho.CP.*12e-6;

% Daily fish production in mg C m^(-2) d^(-1):
FishProduction = 24*sol.Production(end,7).*Parameters.rho.CP.*12e-6;

% Total fish biomass in water column in g C m^(-2):
FishBiomass = trapz(sol.zspan,sol.u(:,7)).*Parameters.rho.CP.*12e-6; 

%% Euphotic and epipelagic dephts

%%% Defining the euphotic depth: 

% Definition of euphotic depth, as fraction of surface light (0.1%,
% Buesseler et.al 2020):
z_euphotic = 0.1e-2*Parameters.E_0;

% The irradiance at each depth in zspan:
E = irradiance(sol.un,Parameters.rho,Parameters.K,Parameters.k, ...
    Parameters.E_0,sol.zspan);

% The euphotic depth is set to the first element in zspan, that is greater
% or equal to the euphotic depth as defined above:
eupthotic_idx = find(E<=z_euphotic,1); % Index at wich E <= z_euphotic.
z_euphotic = sol.zspan(eupthotic_idx);

%%% Defining the epipelagic depth: 

z_epipelagic = 200; % The bottom depth of the epipelagic zone.
% Finding the Index at wich sol.zpsan = z_epipelagic:
[minDiff,epipelagic_idx] = min(abs(sol.zspan-z_epipelagic)); 

% If the value of z_epipelagic is not included in sol.zspan, then a warning
% is issued, stating that the value closest to it will be used instead: 
if minDiff>0
    warning(['The depth %dm is not in the vector zspan, using the ' ...
        'depth closest to %dm (with a difference in value of %gm) ' ...
        'which is %dm.'],z_epipelagic,z_epipelagic, ...
        minDiff,sol.zspan(epipelagic_idx))
end

%% Gravity and diffusive fluxes in mg C m^(-2) h^(-1)

% Using the function ingestionRates to calculate the ingestion rates for H,
% C, Z and F:
I = ingestionRates(sol.un,Parameters.I,Parameters.H,Parameters.sigma);

if optional.NoDVM
    % If there is no diel vertical migration for fish, then the water
    % column ingestion is unaltered:
    expectedIngestion = I.F;
else
    % Calculates the expected value of the ingestion rate for fish in the
    % water column:
    integrand = I.F.*Observations.fishDistribution; 
    expectedIngestion = trapz(sol.zspan,integrand);
end

% Using the function growthLimitation to calculate the limitation fractions
% from equations S16 to S18, and retrieving the logical limitation indices: 
[limFrac,limDepth] = growthLimitation(sol.un,E,Parameters.H, ...
    Parameters.limFrac,Parameters.limDepth);

% Using the function growthRates to calculate the growth rates for
% B, A, and D:
mu = growthRates(Parameters.mu,limFrac,limDepth);

% Using the function Growth to calculate the growth in bacteria, a.
% flagellates, diatoms, h. flagellates, ciliates, mesozooplankton, and
% fish:
growth = Growth(sol.un,mu,Parameters.Y,I,expectedIngestion,N_z);

% Using the function metabolicLossAllocation to calculate the
% allocation of metabolic losses:
[P,R,DOC_HCZ,DOC_F] = metabolicLossAllocation(sol.un,growth, ...
    Parameters.Y,I,expectedIngestion,Parameters.P,Parameters.R, ...
    Parameters.f,Parameters.rho);

% Using the function netGrowth to calculate the biological sources
% and sinks term:
s = netGrowth(sol.un,sol.growth,mu,Parameters.Y, ...
    I,Parameters.f,Parameters.delta,Parameters.k, ...
    P,R,Parameters.rho,DOC_HCZ,DOC_F);

% Extracting kappa and nu from the Parameters structure: 
kappa = Parameters.kappa;
nu = Parameters.nu;

% Using the function shift, to define the matrix containing the elements of
% un shifted up one step, while accounting for equal diffusion from below
% bottom boundary to z = z_max as from z = z_max to below bottom boundary,
% except for P and S, for wich there is a constant supply from below the
% bottom boundary:
[u_shiftUp] = shift(sol.un,s,N_vars,Parameters.P,Parameters.S_b,kappa);

% The net diffusive flux between the the depth z and z + dz:
netDiffusiveFlux = -dz.^(-2).*(kappa.*(u_shiftUp - sol.un))*1e3;

% The downward gravitational flux at depth z:
gravitationalFlux = dz.^(-1).*(nu.*sol.un)*1e3;

% Particle gravity flux, from nmol-P m^(-2) h^(-1) to mg C m^(-2) h^(-1):
particleGravityFlux ...
    = sum(gravitationalFlux(:,8:9),2).*Parameters.rho.CP*12e-6;

% Net diffusive Particulate vertical flux in mg C m^(-2) h^(-1):
particleDiffusiveFlux ...
    = netDiffusiveFlux(:,1).*Parameters.rho.B*12e-6 ...
    + sum(netDiffusiveFlux(:,[2:6 8:10]),2).*Parameters.rho.CP*12e-6;


% Net diffusive DOC vertical flux in mg C m^(-2) h^(-1):
DOCdiffusiveFlux = netDiffusiveFlux(:,12)*12e-6;

%% Active flux in mg C m^(-2) h^(-1)

% Calculating fish processes from the final time step:
fishIngestion = Parameters.rho.CP*12*(I.F.*sol.un(:,7));
fishRespiration = 12*R.F + 12*DOC_F;
fishDefacation = Parameters.rho.CP*12*P.F;
fishMortality = Parameters.rho.CP*12*(Parameters.delta.F*sol.un(:,7));

% The integrands for each depth z=x (0<=x<=z_max) in vector form:
integrands = fishDefacation + fishRespiration ...
    + fishMortality - fishIngestion;

% The active flux for each depth z=x (0<x<=z_max):
activeFlux = -1e-3*flip(cumtrapz(flip(sol.zspan),flip(integrands)));

%% Total flux in mg C m^(-2) h^(-1)
totalFlux = particleGravityFlux + particleDiffusiveFlux ...
    + DOCdiffusiveFlux + activeFlux; % mg C m-2 h-1

%% Relative fluxes (dimensionless)

% The total carbon fluxes: 
totalGravFlux = sum(gravitationalFlux(:,8:9),2) ...
    .*Parameters.rho.CP*12e-6; % Gravity flux in mg C m^(-2) h^(-1). 
totalDiffFlux = particleDiffusiveFlux + DOCdiffusiveFlux; % Diffusive flux.
totalActFlux = activeFlux; % Active flux (redefined only for readability). 

% Defining a logical index vector where all non-zero elements corresponds
% to negative fluxes (i.e., fluxes that are less than zero):
NegDepths = any([totalGravFlux totalDiffFlux totalActFlux] < 0,2);

% Defining a vector of indices, where the first element is the depth index
% following immediately after the last depth index for which any of the
% fluxes are negative (i.e., if "j" is the last depth index for which any
% of the fluxes are negative, then all fluxes have non-negative elements
% for the indices "j+1" to "end"):
nonNegDepths = (find(NegDepths,1,'last')+1):N_z; 

if isempty(nonNegDepths)
    warning(['One or more of the relative fluxes have only negative ' ...
        'values, the relative fluxes are therefore plotted for all ' ...
        'depths in the discretised water column.'])
    nonNegDepths = 1:N_z; 
end

% The relative contribution of gravity flux: 
relativeGravFlux = totalGravFlux(nonNegDepths)./totalFlux(nonNegDepths); 

% The relative contribution of diffusive flux: 
relativeDiffFlux = totalDiffFlux(nonNegDepths)./totalFlux(nonNegDepths); 

% The relative contribution of active flux:
relativeActFlux = totalActFlux(nonNegDepths)./totalFlux(nonNegDepths); 

%% Fluxes at the bottom of the epipelagic zone

% The total carbon export:
TotalCarbonExportEpipelagic = totalFlux(epipelagic_idx)*24;

%  The carbon transfer efficieny (%):
CarbonTransferEfficiencyEpipelagic = TotalCarbonExportEpipelagic/NPP*100;

% The gravitational flux:
GravFluxEpipelagic = particleGravityFlux(epipelagic_idx)*24;

% The diffusive particle flux:
ParticleDiffusiveFluxEpipelagic ...
    = particleDiffusiveFlux(epipelagic_idx)*24;

% The diffusive DOC flux:
DOCdiffusiveFluxEpipelagic = DOCdiffusiveFlux(epipelagic_idx)*24;

% The active flux:
ActiveFluxEpipelagic = activeFlux(epipelagic_idx)*24;

% C:P of export
CPratioExport =  TotalCarbonExportEpipelagic ...
    ./(-24*netDiffusiveFlux(epipelagic_idx,11)*12e-6 ...
    - 24*activeFlux(epipelagic_idx)./(Parameters.rho.CP*1e3));

% The net upward DIP flux in mmol-P m^(-2) d^(-1):
NetUpwardDIPfluxEpipelagic...
    = - netDiffusiveFlux(epipelagic_idx,11)*24e-6;

% The total P export in mmol P m^(-2) d^(d-1):
TotalPexportEpipelagic = -24*netDiffusiveFlux(epipelagic_idx,11)*1e-6 ...
    - 24*activeFlux(epipelagic_idx)./(Parameters.rho.CP*12e3);

%%% The length scale of the gravitational flux attenuation below 200m:
% The length scale of the gravitational flux attenuation below 200m is
% found by solving m equations (where m is the length of
% sol.zspan(epipelagic_idx:end), with two unknowns x1 and x2, on the
% following form: x1 + z(i)x2 = log(particleGravityFlux(i)), where i =
% 1,2,...,m, and x2 is used to calculate the length scale of the
% gravitational fluxt attenuation below 200m.

% The natural logarithm of the particle gravity flux for z in [200,700]:
b = log(particleGravityFlux(epipelagic_idx:end));

% The coefficient matrix A:
if iscolumn(sol.zspan)
    A = [ones(length(b),1) sol.zspan(epipelagic_idx:end)];
else
    A = [ones(length(b),1) sol.zspan(epipelagic_idx:end)'];
end

% Solving the system of linear equations:
x = A\b;

% Calculating the length scale:
LengthScaleGravFlux = -1/x(2);

%% Fluxes at the bottom of the mesopelagic zone

% The gravitational flux:
GravFluxBottom = particleGravityFlux(end)*24;

%  The carbon transfer efficieny:
CarbonTransferEfficiencyBottom = 24*totalFlux(end)/NPP*100;

%% Community respiration
respiration = (R.B + R.HCZ + R.F); % Respiration of B, H, C, Z, and F.

%% Respiration and sequestration of carbon in the mesopelagic

% Integrated mesopelagic respiration in mg C m-2 d-1:
RespirationMesopelagic = 24*(trapz(sol.zspan(epipelagic_idx:end), ...
    respiration(epipelagic_idx:end)*12e-3) + totalFlux(end));

% Weighted integrated mesopelagic respiration: 
wRespiration = sol.zspan.*respiration*12e-3; 
wRespirationMesopelagic = trapz(sol.zspan(epipelagic_idx:end), ...
    24*wRespiration(epipelagic_idx:end)) + 24*700*totalFlux(end);

% Weigthed mean depth of community respiration (m):
WMDR = wRespirationMesopelagic/RespirationMesopelagic; 

% Sequestration time (years) as a function of weigthed mean depth of
% community respiration (m) (regression estimated from Boyd et al 2019,
% fig. 2b):
SequestrationTime = 0.35*WMDR - 12.9;

% The sequestration proxy, in kg C m-2, defined as the product of the
% integrated mesopelagic respiration and the sequestration time:
SequestrationProxy = 365*RespirationMesopelagic*1e-6*SequestrationTime;

%% Input for the time-series of mean concentrations plots

% Collecting plotData in a matrix:
meanConcentrations = [sol.uMean sol.TotalPmean];

% Putting the N_rows by 1 vectors in meanConcentrations in cells and
% storing them in a cell called meanConcentrations:
[N_rows,N_cols] = size(meanConcentrations);
meanConcentrations = cellfun(@(x) mat2cell(x,N_rows,1), ...
    mat2cell(meanConcentrations,N_rows,ones(1,N_cols)), ...
    'UniformOutput',false);

% Defining the integration interval, from t = 0 to t = T:
tspan = days(hours(sol.tspan)); % From hours to days. 

% Storing a copy of tspan in each cell of a cell array of length equal to
% number of tiled plots:
[tSpan{1:N_cols}] = deal({tspan});

% Storing the input in a structure:
inputForPlots.meanConcentrations.tspan = tSpan;
inputForPlots.meanConcentrations.plotData = meanConcentrations;

%% Input for the vertical profiles of solutions plots

% Collecting plotData in a matrix:
verticalProfiles = [sol.u sum(sol.u(:,1:11),2)]; % u and TotalP.

% Putting the N_rows by 1 vectors in verticalProfiles in cells and storing
% them in a cell called verticalProfiles:
[N_rows,N_cols] = size(verticalProfiles);
verticalProfiles = cellfun(@(x) mat2cell(x,N_rows,1), ...
    mat2cell(verticalProfiles,N_rows,ones(1,N_cols)), ...
    'UniformOutput',false);

% Storing a copy of the depth span in each cell of a cell array of length
% equal to number of tiled plots:
[zSpan{1:N_cols}] = deal({sol.zspan});

% Storing input in a structure:
inputForPlots.verticalProfiles.zspan = zSpan;
inputForPlots.verticalProfiles.plotData = verticalProfiles;

%% Input for vertical profiles of carbon flux and fish processes plots

%%% The carbon flux plot:

% Storing the carbon-fluxes in a cell:
carbonFluxes = {particleGravityFlux,particleDiffusiveFlux, ...
    DOCdiffusiveFlux,activeFlux,totalFlux};

% Storing a copy of the depth span in each cell of a cell array of length
% equal to number of line plots in the carbon-fluxes plot:
zSpan = cell(1,length(carbonFluxes));
[zSpan{1:length(carbonFluxes)}] = deal(sol.zspan);

% Storing the depth spans for the carbon-flux plot:
inputForPlots.carbonFluxAndFishProcesses.zspan{1} = zSpan;

% Storing the plotData for the carbon-flux plot:
inputForPlots.carbonFluxAndFishProcesses.plotData{1} = carbonFluxes;

% Values for the horizontal lines in the carbon-fluxes plot:
inputForPlots.carbonFluxAndFishProcesses.z_euphotic = z_euphotic;
inputForPlots.carbonFluxAndFishProcesses.z_epipelagic = z_epipelagic;

%%% The fish processes plot:

% Storing the fish processes in a cell:
fishProcesses = {log10(fishIngestion),log10(fishRespiration), ...
    log10(fishDefacation+fishMortality),log10(respiration*12)};

% Storing a copy of the depth span in each cell of a cell array of length
% equal to number of line plots in the fish processes plot:
zSpan = cell(1,length(fishProcesses));
[zSpan{1:length(fishProcesses)}] = deal(sol.zspan);

% Storing the depth spans for the fish processes plot:
inputForPlots.carbonFluxAndFishProcesses.zspan{2} = zSpan;

% Storing the plotData for the fish processes plot:
inputForPlots.carbonFluxAndFishProcesses.plotData{2} = fishProcesses;

%% Input for vertical profile of the relative carbon flux plot

% Storing the relative carbon-fluxes in a cell:
relativeCarbonFluxes = {relativeGravFlux,relativeDiffFlux,relativeActFlux};

% Storing a copy of the depth span in each cell of a cell array of length
% equal to number of line plots in the carbon-fluxes plot:
zSpan = cell(1,length(relativeCarbonFluxes));
[zSpan{1:length(relativeCarbonFluxes)}] = deal(sol.zspan(nonNegDepths));

% Storing the depth spans for the carbon-flux plot:
inputForPlots.relativeCarbonFluxes.zspan{1} = zSpan;

% Storing the plotData for the carbon-flux plot:
inputForPlots.relativeCarbonFluxes.plotData{1} = relativeCarbonFluxes;

% Values for the horizontal lines in the carbon-fluxes plot:
inputForPlots.relativeCarbonFluxes.z_euphotic = z_euphotic;
inputForPlots.relativeCarbonFluxes.z_epipelagic = z_epipelagic;

%% Input for the observations and simulations plot

%%% The light penetration plot:

% Observed light penetration:
load OBS_den_chl_oxy_pen.mat pen

% Simulated light as the 10-log of fraction of surface light:
Esim = log10(E/Parameters.E_0);

% Storing the depth spans and plot data for the light penetration plot:
inputForPlots.ObservationsAndSimulations.x{1} ...
    = {sol.zspan,sol.zspan(2:end)};
inputForPlots.ObservationsAndSimulations.PlotData{1} = {Esim,pen};

%%% The dissolved inorganic phosphorus (DIP) plot:

% Observed DIP and depths of observations:
load 'OBS_odep,odip,odoc,obac.mat' odep odip

% Dissolved inorganic phosphorus in microM:
DIP = sol.u(:,11)*1e-3;

% Storing the depth spans and plot data for DIP plot:
inputForPlots.ObservationsAndSimulations.x{2} = {sol.zspan,odep};
inputForPlots.ObservationsAndSimulations.PlotData{2} = {DIP,odip};

%%% The silicate (S) plot:

% Observed silicate and observation depths:
odep = [5 10 20 25 40 50 60 76 100 160 200 300 360 400:50:600 690];
oS = [repelem([0.3 0.4 0.5],2) 0.6 0.7 0.8 3 4.8 7.2 ...
    7.7 7.9 8.2 7.2 9.7 10.8 9.6];

% Silicate in microM:
S = sol.u(:,13)*1e-3;

% Storing the depth spans and plot data for the silicate plot:
inputForPlots.ObservationsAndSimulations.x{3} = {sol.zspan,odep};
inputForPlots.ObservationsAndSimulations.PlotData{3} = {S,oS};

%%% The dissolved organic carbon (DOC) plot:

% Observed dissolved organic carbon and observation depths:
load 'OBS_odep,odip,odoc,obac.mat' odep odoc

% DOC in micrroM:
DOC = 50 + sol.u(:,12)*1e-3; % 50 microM is assumed refractory.

% Storing the depth spans and plot data for the DOC plot:
inputForPlots.ObservationsAndSimulations.x{4} = {sol.zspan,odep};
inputForPlots.ObservationsAndSimulations.PlotData{4} = {DOC,odoc};

%%% The dissolved oxygen plot:

% Observed oxygen:
load OBS_den_chl_oxy_pen.mat oxy

% Dissolved oxygen in mL/L:
O = sol.u(:,15)/44.7e3;

% Storing the depth spans and plot data for the oxygen plot:
inputForPlots.ObservationsAndSimulations.x{5} = {sol.zspan,sol.zspan(2:end)};
inputForPlots.ObservationsAndSimulations.PlotData{5} = {O,oxy};

%%% The autotrophs plot:

% Storage matrices:
relChl = nan(N_z,1);
relA = nan(N_z,1);

% The depths at wich the data will be plotted:
plotDepths = 1:181;

% Observed chlorophyll:
load OBS_den_chl_oxy_pen.mat chl
relChl(plotDepths) = ...
    chl(plotDepths)/(trapz(sol.zspan(plotDepths),chl(plotDepths)));

% Autotrophs:
relA(plotDepths) = sum(sol.u(plotDepths,2:3),2) ...
    /(trapz(sol.zspan(plotDepths),sum(sol.u(plotDepths,2:3),2)));

% Storing the depth spans and plot data for the autotrophs plot:
inputForPlots.ObservationsAndSimulations.x{6} = {sol.zspan,sol.zspan};
inputForPlots.ObservationsAndSimulations.PlotData{6} = {relA,relChl};

%%% The bacteria plot:

% Observed bacteria and observation depths:
load 'OBS_odep,odip,odoc,obac.mat' obac odep
relBac = obac/trapz([0; odep; sol.zspan(end)],[obac(1); obac; obac(end)]);

% Bacteria in numbers per mL:
relB = sol.u(:,1)/(trapz(sol.zspan,sol.u(:,1)));

% Storing the depth spans and plot data for the autotrophs plot:
inputForPlots.ObservationsAndSimulations.x{7} = {sol.zspan,odep};
inputForPlots.ObservationsAndSimulations.PlotData{7} = {relB,relBac};

%%% The mesozooplankton plot:

% Observations:

% Depths for zooplankton sampling (Dypvik et al):
zoodepth=[25 75 150 300 500 600];

% Original counts (in m3) (Dypvik et al):
zoo=[601 385 64 35 19 6];

% Relative depth distribution (extrapolated from 25m to surface):
relZoo=zoo/(trapz(zoodepth,zoo)+26*601);

% Mesozooplankton:
relZ = sol.u(:,6)/(trapz(sol.zspan,sol.u(:,6)));

% Storing the depth spans and plot data for the mesozooplankton plot:
inputForPlots.ObservationsAndSimulations.x{8} = {sol.zspan,zoodepth};
inputForPlots.ObservationsAndSimulations.PlotData{8} = {relZ,relZoo};

%%% The gravitational POC flux plot:

% Observations from Torfstein et.al (Summer average, table 3):
gravityFlux = [4.3 1.6 1.3 1.6 1.5];
standardDeviations = [2.4 0.76 0.22 0.44 0.55];
depths = [120 220 350 450 570];

% Storing the depth spans and plot data for the gravitational POC flux
% plot:
inputForPlots.ObservationsAndSimulations.x{9} ...
    = {sol.zspan,sol.zspan,depths};
inputForPlots.ObservationsAndSimulations.PlotData{9} ...
    = {particleGravityFlux,particleGravityFlux+particleDiffusiveFlux,...
    gravityFlux};

% Storing the data for the error bars:
inputForPlots.ObservationsAndSimulations.gravityFlux = gravityFlux;
inputForPlots.ObservationsAndSimulations.depths = depths;
inputForPlots.ObservationsAndSimulations.standardDeviations ...
    = standardDeviations;

%%% The detritus plot:

% Observations from Baker et.al (2017); the amount of non-, slow-, and
% fast-sinking detritus, each amount relative to the total amount:
DetObs50=[0.9493 0.0474 0.00398]; ... at 30-70m depth,
    DetObs150=[0.9407 0.0554 0.00653]; ... at 130-170m depth.
    DetObs = [DetObs50; DetObs150];
Depths = [50; 150]; % Mean observation depths.

% Simulations:
totalDetritus = sum(sol.u(:,8:10),2); % Total detritus in nM-P.
DetSim = sol.u(:,8:10)./totalDetritus; % Relative to total detritus.
DetSim = [DetSim(:,3) DetSim(:,1:2)]; % [non slow fast].

% Storing the depth spans and plot data for the detritus plot:
inputForPlots.ObservationsAndSimulations.x{10} ...
    = {sol.zspan,Depths,sol.zspan,Depths,sol.zspan,Depths};
inputForPlots.ObservationsAndSimulations.PlotData{10} ...
    = {DetSim(:,1),DetObs(:,1),DetSim(:,2),DetObs(:,2), ...
    DetSim(:,3),DetObs(:,3)};

%% Table with additional data:

% Storing the descriptions of the data's ecosystem properties in a
% row-vector for table input:
EcosystemProperty = split("NPP,Euphotic depth," + ...
    "Fish biomass,Fish production,Fluxes at " + string(z_epipelagic) +"m:" ...
    + ",Carbon transfer efficiency," + "Total C export,    Gravitational," + ...
    "    Diffusive (part.)," + "    Diffusive (DOC),    Active," + ...
    "C:P of export,Net upward DIP flux," + "Total P export," + ...
    "Length scale of gravitational flux attenuation below " ...
    + string(z_epipelagic) + "m,Fluxes at " + string(sol.zspan(end)) ...
    + "m:,Carbon transfer efficiency,Gravitational,Respiration and " + ...
    "sequestration of carbon below "  + string(z_epipelagic) + "m:," + ...
    "Respiration,WMDR,Sequestration proxy",',');

% Storing the data units in a row-vector for table input:
Unit = split("mg C m-2 d-1,m,g C m-2," + ...
    "mg C m-2 d-1, ,%," ...
    + strjoin(repelem("mg C m-2 d-1",5),',') + ",mol C/mol P," ...
    + "mmol P m -2 d-1,mmol P m-2 d-1,m, ,%," ...
    + "mg C m-2 d-1, ,mg C m-2 d-1,m,kg C m-2",',');

% Storing the data values in a row-vector for table input:
values = [NPP; z_euphotic; FishBiomass; FishProduction; nan; ...
    CarbonTransferEfficiencyEpipelagic; TotalCarbonExportEpipelagic; ...
    GravFluxEpipelagic; ParticleDiffusiveFluxEpipelagic; ...
    DOCdiffusiveFluxEpipelagic; ActiveFluxEpipelagic; ...
    CPratioExport; NetUpwardDIPfluxEpipelagic; ...
    TotalPexportEpipelagic; LengthScaleGravFlux; nan; ...
    CarbonTransferEfficiencyBottom; GravFluxBottom; nan; ...
    RespirationMesopelagic; WMDR; SequestrationProxy];

% Storing additional data in a table, along with their units and
% description of the ecosystem property the data represents:
SimulationType = sprintf(['kappa=' num2str(kappa(1)/3600,'%g') ...
    '\ndelta.F=%g'],Parameters.delta.F*24*365);
Table = table(EcosystemProperty,Unit,values, ...
    'VariableNames',{'Ecosystem property','Unit',SimulationType});



end

