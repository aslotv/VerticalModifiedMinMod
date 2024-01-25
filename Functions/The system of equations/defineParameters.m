% Defines the parameters in tables S5 and S6.
function [alfa,sigma,mu,Y,I,H,limFrac,limDepth,f,delta,k,nu,K,P,R,rho, ...
    E_0,kappa,S_b] = defineParameters(N_z,options)
%DEFINEPARAMETERS is a function that defines the parameters in the tables
%S5 and S6, in the Supporting Information (Aksnes et al., 2023), results of
%expressions in the system of equations that remains constant in time, and
%storage arrays for results that are not. 
%
%   Required Input: 
%
%   N_z - Number of depths in the discretised water column.
%
%   Optional Input:
%
%   NoFish (NoFish = false as default) - NoFish = true results in a
%   simulation in which fish is not included.
%
%   DiffusivityFactor (DiffusivityFactor = 1 as default) - Will set kappa =
%   DiffusivityFactor*Baseline_kappa.
%
%   FishMortalityFactor (FishMortalityFactor = 1 as default) - Will set
%   delta.F = FishMortalityFactor*Baseline_delta.F.
%
%   TempObsFileName (TempObsFileName = 'temperatureProfile_0to700m.mat' as
%   default) - The name of the MAT-file containing the observed
%   temperatures.
%
%   TempObsVarName (TempObsVarName = 'temp' as default) - The name of the
%   variable containing the observed temperatures which is in the MAT-file
%   with the name given in TempObsVarName.
%
%
%   Output:
%
%   alfa - A structure containing the affinities and clearance rates.
%
%   sigma - Mesozooplankton selectivity factor for C relative to D.
%
%   mu - A structure containing the maximum growth rates for B, A, D, H, C,
%   Z, and F, as well as storage vectors for the growth rates for B, A, and
%   D.
%
%   Y - A structure containing the yields of H, C, Z, F, and B.
%
%   I - A structure containing the maximum ingestion rates of H, C, Z, and
%   F, as well as storage vector for the ingestion rates, and the
%   denominator from the expression defining them.
%
%   H - A structure containing the half-saturation constants for B, A, D,
%   H, C, Z, and F. 
%
%   limFrac - A structure containing storage vectors for the limitation
%   fractions in equations S16 to S19 in the supplementary information. 
%
%   limDepth - A structure containing storage vectors for logical index
%   vectors, where each non-zero elements corresponds to the indices where
%   the corresponding growth limitation is the most limiting. 
%
%   f - A structure containing photosynthetic carbon overflow factor, and
%   the fraction of total loss that goes to detritus and that is respired.
%
%   delta - A structure containing the mortality rates for D and F, and for
%   Z if NoFish = true.
%
%   k - A structure containing the dissolution, leakage and fragmentation
%   rates, as well as empirical coefficients related to chlorophyll light
%   absorption.
%
%   nu - A structure containing the non-zero sinkning velocities.
%
%   K - A structure containing the light attenuation due to clear water and
%   due to other factors than clear water and chlorophyll concentration.
%
%   P - A structure containing the photosynthetic quotient, inorganic
%   phosphate concentration at deep boundary (z_max), and the particulate
%   losses.
%
%   R - A structure containing the respiratory quotient, and respiration
%   losses.
%
%   rho - A structure containing the stoichiometric ratios and conversion
%   factors.
%
%   E_0 - The surface irradiance.
%
%   kappa - The turbulent diffusivity.
%
%   S_b - The inorganic silicate concentration at z = z_max.
%
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%% Arguments block

arguments
    N_z {mustBeInteger,mustBeScalarOrEmpty,mustBeNonempty,mustBePositive}
    options.NoFish {mustBeNumericOrLogical} = false;
    options.DiffusivityFactor (1,1) double = 1;
    options.FishMortalityFactor (1,1) double = 1;
    options.TempObsFileName {mustBeFileName(options.TempObsFileName)} ...
        = 'temperatureProfile_0to700m.mat';
    options.TempObsVarName {mustBeTextScalar} = 'temp';
end

% Extracts variables from structure:
noFish = options.NoFish;
DiffusivityFactor = options.DiffusivityFactor; 
FishMortalityFactor = options.FishMortalityFactor; 
TempObsFileName = options.TempObsFileName;
TempObsVarName = options.TempObsVarName;

%% The maximum growth rates

% Maximum growh rates:
mu.Bmax = 0.25; % Maximum growth rate bacteria.
mu.Amax = 0.054; % Maximum growth rate a. flagellates.
mu.Dmax = 0.06; % Maximum growth rate diatoms.
mu.Hmax = 0.132; % Maximum growth rate h. flagellates.
mu.Cmax = 0.045; % Maximum growth rate ciliates.
mu.Zmax = 0.00625; % Maximum growth rate mesozooplankton.
mu.Fmax = 0.00225; % Maximum growth rate fish.


% The temperature dependies:
mu.T = {'Bmax','Amax','Dmax','Hmax','Cmax','Zmax'};

mu.Q10 = 1.9; % Temperature dependency for maximum growth rate.

%% The affinities and clearance rates
% The variable for the affinities is named alfa, as alpha is the name
% of a function in MATLAB(R).

alfa.BP = 0.08;        % Bacterial affinity for P.
alfa.BL = 1.6*10^(-7); % Bacterial affinity for L.

alfa.AP = 0.04;   % A. flagellate affinity for P.
alfa.AE = mu.Q10*mu.Amax/20;  % A. flagellate affinity for E.

alfa.DP = 0.03;   % Diatom affinity for P.
alfa.DS = 0.001875; % Diatom affinity for S.
alfa.DE = mu.Q10*mu.Dmax/20;  % Diatom affinity for E.

alfa.H = 0.0015;  % H. flagellate clearance rate for B.
alfa.C  = 0.00045; % Ciliate clearance rate for A and H.
alfa.Z  = 0.00015; % Mesozooplankton clearance rate for D.
alfa.F  = 0.00037;  % Fish clearance rate for Z.

sigma = 2; % Mesozooplankton selectivity factor for C relative to D.

% Temperature dependencies:
alfa.T = {'BP','BL','AP','DP','DS','H','C','Z'};

alfa.Q10 = 1.4; % Temperature dependency for nutrient uptake.


%% Temperature dependent affinities and growth rates

%%% Loading the observed temperatures: 
% Loads the observed temperatures and stores them in a structure:
TempObs = load(TempObsFileName,TempObsVarName);

% Extracts the observed temperatures from the structure:
temp = TempObs.(TempObsVarName);

%%% Calculating temperature dependent affinities and maximum growth rates:

% Maximum growth rates with temperature dependency: 
for i=1:length(mu.T) 
    mu.(mu.T{i})=mu.(mu.T{i}).*mu.Q10.^((temp-17)/10);
end

% Affinities with temperature dependency: 
for i=1:length(alfa.T) 
    alfa.(alfa.T{i}) = alfa.(alfa.T{i}).*alfa.Q10.^((temp-17)/10);
end

clear temp TempObs

%% The yields
Y.H  = 0.4;   % H. flagellates yield.
Y.C  = 0.3;   % Ciliates yield.
Y.Z  = 0.15;  % Mesozooplankton yield.
Y.F  = 0.15;  % Fish yield.
Y.BL = 0.001; % Bacterial yield on L.

%% The maximum ingestion rates

% Names of the fields corresponding to h. flagellates, ciliates,
% mesozooplankton, and fish:
fieldNames = ["H", "C", "Z", "F"];

% Maximum ingestion rate for h. flagellates, ciliates, mesozooplankton, and
% fish:
for i = 1:length(fieldNames)
    
    I.(fieldNames(i) + "max") = ...
        mu.(fieldNames(i) + "max")/Y.(fieldNames(i));

    % Note: If fieldNames = ["H", "C", "Z", "F"], then 
    %
    %    I.(fieldNames(1) + "max") = ...
    %        mu.(fieldNames(1) + "max")/Y.(fieldNames(1));
    %
    % is the same as 
    %
    %    I.Hmax = mu.Hmax/Y.H; 

end

%% Detritus sinking and remineralisation
delta.D = 0.01/24; % Diatom mortality (transfer to Det_s and S_opal).

k.opal = 0.25/24; % Dissolution rate of S_opal.
k.l = 0.02/24; % Leakage rate from suspended detritus to P and L.
k.frag = 0.3/24; % Fragmentation rate of Det_s and Det_f.

nu.Det_s = 10/24; % Sinking speed of slow sinking detritus.
nu.Det_f = 100/24; % Sinking speed of fast sinking detritus.

%% Losses to detritus
f.coc = 1; % Photosynthetic carbon overflow.
f.d = 0.2; % Fraction of total loss that enters detritus.

%% Light attenuation coefficients
K.W = 0.00885; % Attenuation from pure water.

% Empirical coefficients related to chlorophyll light absorption:
k.p = 0.10963;
k.e = 0.67175;

%% Oxygen production and consumption
P.Q  = 1; % Photosynthetic quotient.
R.Q  = 1; % Respiratory quotient.

f.rHCZ = 0.6; % Fraction of total loss that is respired (plankton).
f.rF = 0.4; % Fraction of total loss that is respired (fish).

%% Stoichiometric ratios and conversion factors
rho.B = 50; % Molar carbon:phosphorus ratio in B.
rho.CP = 106; % Molar carbon:phosphorus ratio.
rho.DS = 16; % Molar silicate:phosphorus ratio in D.
rho.ChlP = 1/63; % Ratio between Chl and phosphorus in A and D.

%% Forcing variables

% Irradiance: 
E_0 = 600; % Irradiance just below surface.
K.other = 0; % Attenuation other than from pure water and Chl.

% Turbulent diffusivity: 
Baseline_kappa = 3600*3e-4; % Baseline diffusivity.
kappa = DiffusivityFactor*Baseline_kappa; % Scaled diffusivity.

% DIP and silicate concentrations at deep boundary: 
P.b = 1; % Inorganic phosphate concentration at deep boundary (z_max).
S_b = 10.7; % Inorganic silicate concentration at deep boundary (z_max).

% Specific mortality rate for fish:
Baseline_deltaF = 2/(24*365); % Baseline fish mortality.
delta.F = FishMortalityFactor*Baseline_deltaF; % Scaled fish mortality.

% If simulation without fish is desired, define non-zero specific mortality
% rate for mesozooplankton:
delta.Z = noFish/(5*7*24); %  1/(5*7*24) if noFish=true, 0 if noFish=false.

%% Defining constants arising from expressions in the system of PDEs

%%% Half-saturation constants: 

% For bacteria:
H.BP = mu.Bmax./alfa.BP;
H.BL = mu.Bmax./alfa.BL;

% For autotrophic flagellates:
H.AP = mu.Amax./alfa.AP;
H.AE = mu.Amax./alfa.AE;

% For diatoms:
H.DP = mu.Dmax./alfa.DP;
H.DS = mu.Dmax./alfa.DS;
H.DE = mu.Dmax./alfa.DE;

% For heterotrophic flagellates:
H.H = I.Hmax./alfa.H;

% For ciliates:
H.C = I.Cmax./alfa.C;

% For mesozooplankton:
H.Z = I.Zmax./alfa.Z;

% For fish:
H.F = I.Fmax./alfa.F;

%% Storage for Calculations in Time Loop

%%% The growth rates for B, A, and D:
mu.B = zeros(N_z,1);
mu.A = zeros(N_z,1);
mu.D = zeros(N_z,1);

%%% The ingestion rates for H, C, Z, and F:

% For the ingestion rates given in the supporting information:
[I.H,I.C,I.Z,I.F] = deal(zeros(N_z,1)); 

% Additional storage vectors for the ciliate and meesozooplankton
% ingestion rates: 

I.Cdenominator = zeros(N_z,1); % Denominator of C ingestion fraction.
I.CA = zeros(N_z,1); % Ciliate ingestion of a. flagellates. 
I.CH = zeros(N_z,1); % Ciliate ingestion of h. flagellates. 

I.Zdenominator = zeros(N_z,1); % Denominator of Z ingestion fraction.
I.ZD = zeros(N_z,1); % Mesozooplankton ingestion of diatoms.  
I.ZC = zeros(N_z,1); % Mesozooplankton ingestion of ciliates. 


%%% The nutrient and light limitation fractions:

% For bacteria: 
limFrac.BP = zeros(N_z,1); 
limFrac.BL = zeros(N_z,1); 

% For autotrophic flagellates: 
limFrac.AP = zeros(N_z,1); 
limFrac.AE = zeros(N_z,1);

% For diatoms: 
limFrac.DP = zeros(N_z,1);
limFrac.DS = zeros(N_z,1); 
limFrac.DE = zeros(N_z,1); 


%%% Logical nutrient and light limitation depth indices:

% For bacteria: 
limDepth.BP = false(N_z,1); % Depths where DIP is most limiting.
limDepth.BL = false(N_z,1); % Depths where l-DOC is most limiting. 

% For autotrophic flagellates: 
limDepth.AP = false(N_z,1); % Depths where DIP is most limiting. 
limDepth.AE = false(N_z,1); % Depths where irradiance is most limiting. 

% For diatoms: 
limDepth.DP = false(N_z,1); % Depths where DIP is most limiting.
limDepth.DS = false(N_z,1); % Depths where silicate is most limiting. 
limDepth.DE = false(N_z,1); % Depths where irradiance is most limiting. 

%%% Particulate losses to detritus:
P.HC = zeros(N_z,1); % From H and C to Det_s.
P.Z  = zeros(N_z,1); % From Z to Det_f.
P.F  = zeros(N_z,1); % From F to Det_f.

%%% Respiration losses:
R.B   = zeros(N_z,1); % From B.
R.HCZ = zeros(N_z,1); % From H, C, and Z.
R.F   = zeros(N_z,1); % From F.
end

function mustBeFileName(FileName)
% This function simply corrects the slashes (\ or /) in the file name given
% as input, if necessary, and passes the file name to the MATLAB(R)
% argument validation function mustBeFileName. If mustBeFile throws an
% error, the function checks if the file is in a folder added to the
% MATLAB(R) search path; if it is, then the argument is valid, if not, the
% error produced by mustBeFile is thrown. 
try
    FileName = correctSlashes(FileName);
    mustBeFile(FileName)
catch ME
    if isempty(which(FileName))
    throwAsCaller(ME);
    end
end
end
