%% Vertical Modified Minimum Model
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

% Note: Enabling code folding is advised, as this may increase the
% readability of the code.

%%% The Vertical Modified Minimum Model
% The 0D minimum microbial food web model (or MinMod) of
% Thingstad et al. (2007,2020), is here set up for a 1D water column.
% Light, zooplanktivorous fish, and three detritus groups have been added
% to the original MinMod.

%% Pre-run preparations

% Gives the option to save workspace as a MAT-file if not empty, and adds
% the folders with the names given as input to the MATLAB(R) search path if
% they're not already on it (note that if any variables are defined in this
% script, above the following line, these will also be saved if any
% pre-existing variables are to be saved):
FolderName = preRunPreparations(who,'Functions','Observations', ...
    'SimulationResults\FinalResults');
clearvars -except FolderName % Clears workspace.

% Defining a string which will be used to display the time at which the
% programme started in the waitbar:
start = sprintf( "Started current run at " ...
    + string(datetime('now','format','h.mma')) + ".");
tic % Starts timing the run.

%% The solution matrix u
% The solutions of the system of partial differential equations, are stored
% in a matrix u, u=[u1 u2 ... u15], where u1,u2,...,u15 are vectors storing
% the solutions of the 15 equations.
%
% Each vector uj has the dimensions (N_z,1), where N_z is the number of
% depth cells.
%
% The solution matrix u has the dimension (N_z,N_vars), where N_vars is
% the number of state variables, or equations, in the system of partial
% differential equations; uj = u(:,j).
%
%%% Correspondence between the vectors uj and the state variables:
% * u1  - B, bacteria,
% * u2  - A, autotrophic flagellates,
% * u3  - D, diatoms,
% * u4  - H, heterotrophic flagellates,
% * u5  - C, ciliates,
% * u6  - Z, mesozooplankton,
% * u7  - F, fish,
% * u8  - Det_s, slow sinking detritus,
% * u9  - Det_f, fast sinking detritus,
% * u10 - Det_n, non-sinking detritus,
% * u11 - P, inorganic phosphate,
% * u12 - L, Labile dissolved organic carbon (L-DOC),
% * u13 - S, silicate,
% * u14 - S_opal, particulate organic silicate,
% * u15 - O, dissolved oxygen.

% Number of variables/equations in the system of PDEs:
N_vars = 15;

%% Spatial and Temporal Discretisation
% The following sections defines the spatial and temporal discretisation of
% the system of PDEs.

%%% Spatial discretisation:
% The depth is given in meters (m).
z_0 = 0; % The minimum depth of the water column*.
z_max = 700; % The maximum depth of the water column.
dz = 1; % Length of depth step**.

% Defining the number of depths (N_z), a vector containing the depths
% (zspan), and the depth range (zrange) in the discretised water column, as
% well as checking the feasibility of the value of the spatial step length
% with respect to the depth interval (z_0,z_max) and returning an alternate
% dz if neccessary (Note that the N_z returned corresponds to the altered
% dz if it is changed):
[N_z,dz,zspan,zrange] = spatialDiscretisation(z_0,z_max,dz);

% * This code is written for z_0 = 0, due to the surface integrals.
% ** All default observations are defined for dz = 1, in order to use a
% different depth step length it is necessary to use observations defined
% for the new depth step length.

%%% Temporal discretisation:
% The time t, is given in hours (h).
t = 0; % Initial time.
T_days = 75000; % The final time of the integration interval (d)
T = hours(days(T_days)); % The final time of the integration interval (h).
dt0 = 1; % Initial temporal step length.

% The t-values for which the results will be stored: 
tspan = hours(days(0:T_days)); % t-values corresponding to the solutions.

% The embedded Runge-Kutta tableu which will be used in the
% function eERK.m:
[c,A,b,d,p] = ERK435Tableu; % Defined by the function "ERK345Tableu.m".

% Checking if the Runge-Kutta method has the first same as last
% (FSAL) property:
FSAL = all(A(end,:) == b) && b(end) == 0;

% Reshaping the A, b, and d from the RK-tableu, to enable
% vectorised code in the embedded RK-functions:
A = reshape(A,length(A(:,1)),1,length(A(1,:)));
b = reshape(b,1,1,length(b));
d = reshape(d,1,1,length(d));

%% Defining the Simulation Scenario
% This section defines variables that ensures the neccessary alterations of
% the code for running the simulated scenarios presented in the manuscript.
% Note that the function "SetScenario.m" also gives the option to set the
% values of these four variables manually (see "SetScenario.m" for
% details). In addition, the names of the MAT-files containing the
% observation variables, as well as the name of the variables, are defined.

% Setting values to variables defining the simulation scenario, and storing
% the chosen scenario:
[NoFish,NoDVM,DiffusivityFactor,FishMortalityFactor,Scenario] = ...
    SetScenario("Baseline");

% In the case of "SetManually" being given as input to SetScenario, the
% user can cancel the execution of SetScenario.m, in this case the
% simulation should be stopped:
if any(cellfun(@isempty,{NoFish,NoDVM,DiffusivityFactor, ...
        FishMortalityFactor,Scenario}))
    warndlg(['Execution of SetScenario.m was canceled and the run ' ...
        'was aborted.'])
    return
end

%%% The observations:
% Note that these observations are not used in the figure where the results
% from the simulation is plotted against what is there reffered to as
% observations; in these simulations/observations plots, the unaltered
% observations are used (see the script "DefiningObservationMATfiles.m"
% and the function "postProcessing.m" for details).

% The temperature profile:
TempObsFileName = 'temperatureProfile_0to700m.mat';
TempObsVarName = 'temp';

% The oxygen profile:
OxyObsFileName = 'oxygenProfile_0to700m.mat';
OxyObsVarName = 'oxy';

% The fish distribution:
FishDistFileName = 'fishdistribution_0to700m.mat';
FishDistVarName = 'fishDistribution';

%%% Switching between expected F-ingestion and depth specific F-ingestion:
% Defining the function handle for calculating fish ingestion and the fish
% distribution variabel based on the value of NoDVM:
if NoDVM
    FishIngestion = @FishIngestionNoDVM;
    fishDistribution = []; 
else
    FishIngestion = @FishIngestion;
    % Defining the variable fishDistribution as it is used in the time
    % loop:
    ObsFishDist = load(FishDistFileName,FishDistVarName);
    fishDistribution = ObsFishDist.(FishDistVarName);
end

%% The Parameters and Initial Conditions
% These sections defines the parameters and the initial conditions, given
% the simulation scenario set up in the previous section.

%%% The Parameters:
% The parameters from tables S5 and S6, and equations S19 to S22 in the
% supporting information (see defineParameters.m for details):
[alfa,sigma,mu,Y,I,H,limFrac,limDepth,f,delta,k,nu,K,P,R,rho, ...
    E_0,kappa,S_b] = defineParameters(N_z, ...
    "NoFish",NoFish, ...
    "DiffusivityFactor",DiffusivityFactor, ....
    "FishMortalityFactor",FishMortalityFactor, ...
    "TempObsFileName",TempObsFileName, ...
    "TempObsVarName",TempObsVarName);

% Putting kappa and nu in vector forms to simplify the diffusive and
% advective terms (i.e., vectorising):
[kappa,nu] = defineDiffSink(kappa,nu,N_vars);

%%% The initial conditions:
% Defining the initial conditions for the current simulation, and storing
% the initial conditions from the initial simulation (i.e., the initial
% conditions that where defined explicitly within the function): (See the
% function "initialConditions.m" for details)
[u,u0,growth0,u0Mean,TotalP0mean,Production0] = initialConditions("set",...
    "N_z",N_z, ...
    "N_vars",N_vars, ...
    "zspan",zspan, ...
    "zrange",zrange, ...
    "FishDistFileName",FishDistFileName, ...
    "FishDistVarName",FishDistVarName, ...
    "OxyObsFileName",OxyObsFileName, ...
    "OxyObsVarName",OxyObsVarName, ...
    "NoFish",NoFish);

%% Storage Matrices for Additional Results from the Time Loop

%%% Temporary storage matrices: 
% The following matrices will store the results for each time step, until
% the number of iterations surpasses the number of hours in the integration
% interval, or the pre-defined maximum number of time steps that will be
% stored temporarily. When this happens, the results are interpolated and
% stored for every 24h, the iteration counter is reset, and the temporary
% storage matrices will be reused.

% Defining the number of time steps that will be stored temporarily as the
% minimum of two times the number of days in the integration interval
% (T_days) and the maximum number of time steps that will be stored
% temporarily:
N_temp = min(2*T_days,20000); 

% The solutions found at each time step: 
u_temp = cat(1,reshape(u,1,N_z,N_vars),zeros(N_temp,N_z,N_vars)); 

% The osmo- and phagotroph growth: 
growth_temp = cat(1,reshape(growth0,1,N_z,7),zeros(N_temp,N_z,7)); 

% The t-values corresponding to the solutions:
tspan_temp = [t; zeros(N_temp,1)]; 

% Solutions averaged over the water column:
uMean_temp = [u0Mean; zeros(N_temp,N_vars)]; 

% Total P averaged over the water column:
TotalPmean_temp = [TotalP0mean; zeros(N_temp,1)];

% Production of B, A, D, H, C, Z and F: 
Production_temp = [Production0; zeros(N_temp,7)]; 

%%% Final storage matrices: 
% The following matrices will store the interpolated values of the
% solutions found while iterating, and correspond to the solutions
% evaluated each time t in tspan (e.g., tspan = [0 24 48 72 96 ...]).
N_t = length(tspan); % The number of t-values for which results are stored.
uMean = zeros(N_t,N_vars); % Solutions averaged over the water column.
TotalPmean = zeros(N_t,1); % Total P averaged over the water column.
Production = zeros(N_t,7); % Production of B, A, D, H, C, Z and F.

%% Preparation for Dynamic Time-Step-Loop
% These sections defines the variables needed for the dynamic
% time-step-loop.

%%% Error control and non-negative solution indices:

Tol = 5e-8; % Error tolerance for local error.

% The indices of the variables to impose non-negative value condition on:
nonNegIndices = 1:14; % All except dissolved oxygen (O).

% Defining the error as a function of the local error and the "negative"
% error (i.e., the measure of how negative the non-negative solutions are):
eFun = @(err,uNonNeg) max(norm(err,'inf'),NegNorm(uNonNeg,'inf'));

% Defining the factor adjusting the length of the time steps, as a function
% of the error and the tolerance (Butcher 2016):
ButcherFactor = @(Tol,e) max(0.5,min(2.0,0.9*nthroot(Tol/e,p+1)));

%%% Defining the variables used to determine if steady state is reached:

% The indices of the variables for which steady-state solutions are
% desired:
steadyIndices = 1:14; % All except dissolved oxygen (O).

% Solutions varying within +- steadyTol for one day of time steps (24
% hours) are accepted as steady state solutions:
steadyTime = 24; % One day in hours (1d=24h).

% The number of decimals to remain constant throughout steadyTime:
steadyDecimals = 4;

% Steady tolerance per time step as a function of dt_n (i.e., the time
% steps):
steadyTol = @(dt_n) (0.5*10^-steadyDecimals)/steadyTime*dt_n;

%%% Defining initial values for counters and logicals:

n = 0; % Counter for number of successful time steps (outer iterations).

d1 = 1; % The initial value for the first index in the quiery points. 

d2 = 2; % The initial value for the last index in the quiery points. 

% Defining a storage variable (or initial value) for the sum of consecutive
% time steps with solutions meeting the steady tolerance:
j = 0; 

% The original value of the final time (t) of the integration interval (T)
% is stored, in case the solutions cease to meet the steady condition while
% finishing integration:
Toriginal = T;

% Defining a minimum threshold for the time step based on the threshold
% in the solver "ode45.m" by MATLAB(R) (type "open ode45" in the
% command window to see the code) and the documentation for spacing of
% floating point numbers (see documentation for eps for details):
dtThresh = @(t) eps(t);

% Defining initial value for the failsafe logical abortIntegration:
abortIntegration = false;

% Initial values for the logicals resulting from the "per-time-step" steady
% check:
steady = false(size(u(:,steadyIndices))); 

%% The Time Loop
% The time loop will iterate as long as t is less than the final
% time T (note that T will be redefined if the solutions are within the
% steady tolerance and the redefined final time T will be less than,
% possibly equal to, the original final time).

% Initialising the progress dialog window (change argument 'on' to 'off',
% to disable the progress dialog window):
timeLoopProgress(start,t,T,'on'); % See timeLoopProgress.m for details. 

% Note that disabling the progress dialog window results in loosing the
% opportunity to (smoothly) stopping the time loop mid integration and
% still saving the results found up until the simulation was canceled.

while t < T

    if exist('w','var') % If the progress dialog window is enabled.
        % Updates the uiprogressdlg:
        w.Value = max(t/T,nnz(steady)/numel(steady)); 

        % Stops while loop if run is canceled:
        if w.CancelRequested; delete(uifig); break; end
    end

    % The counter for successful time steps is increased by 1 (note that
    % this counter is reset if the temporal integration is aborted or
    % canceled):
    n = n+1; 
    
    % The solution matrix from the previous time step is stored in the
    % variable un:
    un = u; 

    % Defines a logical scalar, which remains false until a time step is
    % deemed successful:
    successfulTimeStep = false;

    % The inner loop, call it the error optimisation loop, will run as long
    % as the local error is larger than the local error tolerance, or as
    % long as any of the solutions not allowed to be negative are negative.
    while ~successfulTimeStep

        %%% Making an attempted time step:

        % The initial (or refined) length of the current time step:
        dt_n = dt0; 

        % Makes a time step with the length of the time step defined above:
        [u,err,growth] = eERK(un,dt_n,c,A,b,d,FSAL, ...
            sigma,mu,Y,I,H,limFrac,limDepth,f,delta,k,nu,K,P,R,rho,E_0, ...
            kappa,S_b,fishDistribution,FishIngestion,dz,zspan,N_z,N_vars);

        %%% Calculating the error:

        % Defining a logical column vector, where each non-zero value
        % corresponds to the rows in u(:,nonNegIndices) where one or more
        % of the elements have a negative value:
        negDepths = any(u(:,nonNegIndices)<0,2);
        
        % The error is defined as the maximum of the local error and the
        % maximum value of the absolute value of the negative value of u;
        % the local error is defined as the maximum absolute row sum of the
        % difference between the two solution matrices u and u_t:
        e = eFun(err,u(negDepths,nonNegIndices));

        % The following tolerance is equal to Tol if the local error is
        % larger than Tol or none of the non-negative variables have
        % negative values at any depth; it is equal to zero otherwise:
        nonNegTol = Tol*(e>=Tol || ~any(negDepths));

        % Using a safety factor from Buthcer (2016) to adjust time step:
        dt0 = ButcherFactor(nonNegTol,e)*dt_n;

        % Stopping the inner while loop, and setting the value of
        % abortIntegration equal to true if the adjusted time step is
        % smaller than the threshold:
        if dt0 < dtThresh(t); abortIntegration = true; break; end

        %%% Determining if the time step was successful:

        % Defining a logical scalar, which is true only if the time step
        % was successful:
        successfulTimeStep = e<Tol && ~any(negDepths);
    end

    % Stops the outer while loop if the length of the time steps has fallen
    % below the threshold:
    if abortIntegration; n = n-1; break; end

    % Checks if all elements of the matrix abs(u-un) are less than
    % steadyTol(n):
    steady = abs(u(:,steadyIndices)-un(:,steadyIndices)) < steadyTol(dt_n);
    isSteady = all(steady,'all'); % True if all are steady.

    %%% Storing desired model output:

    if t + dt_n < tspan(d2)

        % The solutions averaged over the entire water column:
        uMean_temp(n+1,:) = trapz(zspan,u)/zrange;

        % Average Total P in water column:
        TotalPmean_temp(n+1) = trapz(zspan,sum(u(:,1:11),2))/zrange;

        % The production is defined as the total growth for each osmo- and
        % phagotroph (in nmol-P m^-2 h^-1):
        Production_temp(n+1,:) = trapz(zspan,growth(:,:))*1e3;

        % Defining the time t(n) by adding the length of the time step just
        % taken, dt_n, to the time t (where t = t(n-1)):
        tspan_temp(n+1) = t+dt_n;

        % Storing the solutions found in the current time step:
        u_temp(n+1,:,:) = u;

        % Storing the osmo- and phagotroph growth found in the current time
        % step:
        growth_temp(n+1,:,:) = growth;

        % Resetting the iteration counter: 
        n = n-1; 
    else

        % Advancing the iteration counter: 
        n = n+1; 

        % The solutions averaged over the entire water column:
        uMean_temp(n+1,:) = trapz(zspan,u)/zrange;

        % Average Total P in water column:
        TotalPmean_temp(n+1) = trapz(zspan,sum(u(:,1:11),2))/zrange;

        % The production is defined as the total growth for each osmo- and
        % phagotroph (in nmol-P m^-2 h^-1):
        Production_temp(n+1,:) = trapz(zspan,growth(:,:))*1e3;

        % Defining the time t(n) by adding the length of the time step just
        % taken, dt_n, to the time t (where t = t(n-1)):
        tspan_temp(n+1) = t+dt_n;

        % Storing the solutions found in the current time step:
        u_temp(n+1,:,:) = u;

        % Storing the osmo- and phagotroph growth found in the current time
        % step:
        growth_temp(n+1,:,:) = growth;

        % Advancing to next quiery point: 
        d2 = d2+1; 
    end

    % If all elements of the matrix abs(u-un) are less or equal to
    % steadyTol, the current time step is added to j; the value of j is
    % reset otherwise:
    j = (j+dt_n)*isSteady; % j = j + dt if isSteady = true, j = 0 if not.

    % Storing the time at wich the steady condition was met (if met):
    if j < steadyTime

        % While j is less than steadyTime, the final time T remains
        % unchanged (i.e., the final time T is as defined in the temporal
        % discretisation):
        T = Toriginal;

        % Clearing the variable storing the time at which the steady
        % condition was met, if the solutions no longer meet said
        % condition:
        clear tSteady

    elseif ~exist("tSteady","var")

        % When the solutions becomes steady (i.e., when j is greater or
        % equal to steadyTime and if the time at which the became steady is
        % not already stored), the final time T is redefined so that the
        % integration continues throughout the day started:
        T = ceil(tspan_temp(n+1)/24)*24;

        % Storing the time at which the steady condition was met:
        tSteady = tspan_temp(n+1);

    end

    %%% Preparing for the next time step:

    % Defining the value of t for the next time step (i.e., in the next
    % iteration of the outer loop, this t is equal to t(n-1) until updated
    % again as below):
    t = t + dt_n; 

    % Adjusting initial value for next time step if it is too big or if
    % last time step will become too small:
    if t + dt0 > T || T - ( t + dt0 ) < eps(T)
        dt0 = T-t;
    end

    %%% Permanently storing results if necessary: 
     if n >= N_temp

         % Resetting the value of d2 if necessary: 
         if tspan_temp(end) < tspan(d2)
             d2 = d2-1; 
         end

         % Defining the quiery points: 
         tq = tspan(d1:d2); 

         % Finding the index of first element that is greater or equal to
         % the last element in the quiery points:
         n_d2 = find(tspan_temp >= tq(end),1,"first"); 

         % Using the MATLAB(R) function interp1, to interpolate the average
         % water column concentrations and the production, and define them
         % as the values at each of the times in tspan:
         uMean(d1:d2,:) ...
             = interp1(tspan_temp(1:n_d2),uMean_temp(1:n_d2,:),tq);
         TotalPmean(d1:d2,:) ...
             = interp1(tspan_temp(1:n_d2),TotalPmean_temp(1:n_d2,:),tq);
         Production(d1:d2,:) ...
             = interp1(tspan_temp(1:n_d2),Production_temp(1:n_d2,:),tq); 

         % Updating the iteration counter and the temporary storage
         % matrices:
         n = length(n_d2:N_temp+1)-1;
         tspan_temp(1:n+1,:) = tspan_temp(n_d2:end); 
         uMean_temp(1:n+1,:) = uMean_temp(n_d2:end,:);
         TotalPmean_temp(1:n+1,:) = TotalPmean_temp(n_d2:end,:);
         u_temp(1:n+1,:,:) = u_temp(n_d2:end,:,:); 
         growth_temp(1:n+1,:,:) = growth_temp(n_d2:end,:,:); 

         % Updating the first index of the quiery points: 
         d1 = d2+1; 

         % Updating the initial value of d2: 
         d2 = d1 + 1; 
    end
end

%%% Interpolating the remaining results (if any): 

% The index of the last element in tspan that is less or equal to the
% current value of t (i.e., tspan_temp(n+1)):
d2 = find(tspan <= tspan_temp(n+1),1,'last');

% Defining the quiery points:
tq = tspan(d1:d2);

if ~isempty(tq)
    % Finding the index of first element that is greater or equal to the
    % last element in the quiery points:
    n_d2 = find(tspan_temp >= tq(end),1,"first");

    % Using the MATLAB(R) function interp1, to interpolate the average
    % water column concentrations and the production, and define them as
    % the values at each whole day:
    uMean(d1:d2,:) ...
        = interp1(tspan_temp(1:n_d2),uMean_temp(1:n_d2,:),tq);
    TotalPmean(d1:d2,:) ...
        = interp1(tspan_temp(1:n_d2),TotalPmean_temp(1:n_d2,:),tq);
    Production(d1:d2,:) ...
        = interp1(tspan_temp(1:n_d2),Production_temp(1:n_d2,:),tq);
end

% Defining the solution matrices u and un, and the growth matrix growth, in
% the event that temoral integration was canceled or aborted:
if ~isempty(tq) && tspan_temp(n+1)~=tq(end)

    % If the temporal integration was canceled or aborted (i.e., if
    % tspan_n(n+1) does not equal tq(end)), then the solution matrices
    % u and un, and growth matrix are defined by interpolation of their
    % respective temporary storage matrices:
    un = ...
        squeeze(interp1(tspan_temp(1:n_d2),u_temp(1:n_d2,:,:),tq(end-1)));
    u = ...
        squeeze(interp1(tspan_temp(1:n_d2),u_temp(1:n_d2,:,:),tq(end)));
    growth = squeeze(interp1(tspan_temp(1:n_d2), ...
        growth_temp(1:n_d2,:,:),tq(end)));
end

% Deletes the uiprogressdlg if it exists:
if exist('uifig','var'); delete(uifig); end

time = toc; % Stops the timer and stores the run time.

%% Post-Processing

%%% Storing information about stop-condition and playing an alert sound:

% Defines tSteady as an empty array if steady state was not reached: 
if ~exist("tSteady","var"); tSteady = []; end

% Using the function determineTerminationCause.m to determnine and store
% the termination cause for the temporal integration, play a sound
% alerting to the fact that the temporal integration is terminated, and
% give the option to save the results of the simulation: 
[stopConditionMet, saveResults] ...
    = determineTerminationCause(j,steadyTime,steadyTol,tSteady,T,t,dt_n,...
    abortIntegration,dtThresh,nonNegTol,e);

%%% Putting the parameters in a structure for saving:

fieldNames = split(['alfa,sigma,mu,Y,I,H,limFrac,limDepth,f,' ...
    'delta,k,nu,K,P,R,rho,E_0,kappa,S_b'],',')';

Parameters = cell2struct({alfa,sigma,mu,Y,I,H,limFrac,limDepth, ...
    f,delta,k,nu,K,P,R,rho,E_0,kappa,S_b},fieldNames,2);

%%% Putting the observation variables in a structure for saving:

% The diel vertical fish distribution:
Observations.fishDistribution = fishDistribution;

% The oxygen profile:
OxyObs = load(OxyObsFileName,OxyObsVarName);
Observations.oxygenProfile = OxyObs.(OxyObsVarName);

% The temperature profile:
TempObs = load(TempObsFileName,TempObsVarName);
Observations.temperatureProfile = TempObs.(TempObsVarName);

%%% Saving, post processing, and plotting results:

% Defining the prefix for the MAT-file:
Version = Scenario; 

% If the solutions met the steady state condition, then the logical Steady
% equals true; if not Steady = false:
Steady = j >= steadyTime; 

% Using the function saveAndPlotResults.m to save (if desired), post
% process and plot the results:
[sol,totalRunTime,RunTime,inputForPlots,Table] = ...
    saveAndPlotResults(initialConditionsInput,Parameters, ...
    Observations,NoDVM,z_0,z_max,dz,u0,un,u,Tol, ...
    uMean(1:d2,:),TotalPmean(1:d2),growth, ...
    Production(1:d2,:),tspan(1:d2),time, ...
    stopConditionMet,j,Version,Steady,Scenario,saveResults);

% Displays the run time in hours, minutes and seconds:
h = time/(3600);
m = (h-floor(h))*60;
s = (m-floor(m))*60;
fprintf(['\nFinding the solutions to the equations took %d hours,' ...
    ' %d minutes and %.3f seconds.\n'], floor(h),floor(m),s)

%% Post-run Clean Up

%%% Cleaning up the workspace: 

% Storing the names of the variables to keep in the workspace:
theseVars = cellstr(split("FolderName,stopConditionMet," + ...
    "initialConditionsInput,Parameters,Observations,sol," + ...
    "totalRunTime,RunTime,inputForPlots,Table,Scenario,Tol",","));

% Clearing all variables in the workspace, except the variables whos names
% are stored in theseVars (i.e., clear all variables except these
% variables):
clearvars('-except',theseVars{:})

%%% Restoring MATLAB(R) search path: 

% Giving the option to remove the folders added to the MATLAB(R) search
% path at the beginning of the programme from the search path again: 
restorePath(FolderName);

% clears the variable FolderName: 
clear FolderName