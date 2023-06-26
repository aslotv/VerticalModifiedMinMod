% This function plots the results. 
function [figures,tilingLayouts,ax] = plotResults(input)
%PLOTRESULTS is a function that plots the results shown in the both the
%publication and the supporting information. 
% 
%   Optional input:
%
%   FileName (FileName = Baseline_steady.mat" as default) - The name of the
%   file in which the results to be plotted are stored.
%
%   InputForPlots - A structure where the input for the plots are stored
%   (given as output argument from the function postProcessing.m).
%
%   Output:
%
%   figures - A structure containing the handels of the figures. 
%
%   tilingLayouts - A structure containing the handels of the tiled
%   layouts. 
%
%   ax - A structure containing the axes of the tiled plots. 
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

arguments
    input.FileName (1,:) {mustBeText} = "Baseline_steady.mat"
    input.InputForPlots {mustBeA(input.InputForPlots,'struct')}
end
%% Defining the Plot Input

if ~isfield(input,'InputForPlots')
    % If the input for the plots was passed directly as input to this
    % function, then the solutions from t=0, to t = T, and the parameters
    % are loaded from the MAT-file with the name stored in FileName:
    load(input.FileName,'inputForPlots')
else 
    % If the input for the plots was passed directly as input to this
    % function, then the variable inputForPlots is defined as the structure
    % given as function input:
    inputForPlots = input.InputForPlots; 
end

%% Color variables as RGB triplets 

% MATLAB(R) default colors for plots (as RGB triplets):
colors = colororder('default'); close(gcf)
plotColors = mat2cell(colors,ones(1,7)); 
[blue,orange,yellow,purple,green,lightBlue,red] = deal(plotColors{:}); %#ok<ASGLU> 

% Standard line color, equivalen to 'k' and 'black':
black = zeros(1,3); 

%% Set up for Figures 1 and 2

% Titles for the tiled plots:
noConversionTitles ...
    = { 'Bacteria (nM-P)', 'A-Flag (nM-P)', 'Diatoms (nM-P)', ...
        'H-Flag (nM-P)','Ciliates (nM-P)', 'Meso-Zooplankton (nM-P)', ...
        'Fish (nM-P)', 'Detritus slow (nM-P)', 'Detritus fast (nM-P)', ...
        'Suspended detritus (nM-P)', 'DIP (nM-P)', 'L-DOC (nM-C)', ...
        'Silicate (nM-Si)', 'POS (nM-Si)', 'Oxygen (nM-O_2)', ...
        'Tot P (nM-P)' ...
        };

% A vector of indices, which will give the desired ordering of the plots:
desiredLayout = [4:6; ... H, C, and Z.
                1:3; ... B, A, and D.
                12 11 13; ... L, P, and S.
                8 15 9; ... Det_s, O, and Det_f.
                16 10 7; ... TotalP, Det_n, and F.
                ];

[m,n] = size(desiredLayout); % Rows and columns for tiling.
%% Figure 1: Time-series of mean concentrations:

tSpan = inputForPlots.meanConcentrations.tspan; 
meanConcentrations = inputForPlots.meanConcentrations.plotData; 

% Using the function plotInTiles to make the tiled plot:
[figures.meanConcentrations, ...
 tilingLayouts.meanConcentrations, ...
 ax.meanConcentrations] ...
    = plotInTiles(tSpan,meanConcentrations, ...
                  'Layout', {m,n}, ...
                  'FigureName', 'MeanConcentrations', ...
                  'PlotOrder', reshape(desiredLayout',1,m*n), ...
                  'Titles',noConversionTitles,...
                  'MainXlabel','Time (in days)', ...
                  'MainYlabel','Concentrations (nM)', ...
                  'MainTitle','Average Water Column Concentrations');

%% Figure 2: Vertical profiles of solutions

zSpan = inputForPlots.verticalProfiles.zspan; 
verticalProfiles = inputForPlots.verticalProfiles.plotData; 

[figures.verticalProfiles, ...
 tilingLayouts.verticalProfiles, ...
 ax.verticalProfiles] ...
= plotInTiles(verticalProfiles,zSpan, ...
              'Layout', {n,m}, ...
              'FigureName', 'VerticalProfiles', ...
              'PlotOrder', reshape(desiredLayout,1,m*n), ...
              'Titles',noConversionTitles, ...
              'MainXlabel',"Concentrations (nM)", ...
              'MainYlabel',"Depth (in meters)", ...
              'MainTitle',"Vertical profiles of steady-state solutions", ...
              'DepthPlot',true);

%% Set up for figure 3

% Defining options for the Vertical Carbon Flux plot:
Legends{1} = {'Gravity POC flux', 'Net diffusive POC flux', ...
              'Net diffusive DOC flux', 'Net active flux', 'Total flux'};
Colors{1} = [red; yellow; green; blue; black]; 
LineSpecs{1} = {'-', '-', '-','-','-'};
Titles{1} = 'Vertical Carbon Flux';
XLabels{1} = 'Net downward flux (mg-C m^{-2} h^{-1})';
YLabels{1} = 'Depth (m)';

% Defining options for the Fish Processes plot:
Legends{2} = {'Fish ingestion', 'Fish excretion and respiration', ...
                  'Fish defacation and mortality', ...
                  'Total community respiration'};
Colors{2} = [blue; green; red; black]; % 4 by 3 matrix of RGB triplets
LineSpecs{2} = {'-', '-', '-','--'};
Titles{2} = 'Fish processes';
XLabels{2} = 'Carbon intake/loss (log_{10}, \mu g-C m^{-3} h^{-1})';
YLabels{2} = 'Depth (m)';
%% Figure 3: Vertical profiles of carbon flux and fish processes

plotData = inputForPlots.carbonFluxAndFishProcesses.plotData; 
zSpan = inputForPlots.carbonFluxAndFishProcesses.zspan; 
z_euphotic = inputForPlots.carbonFluxAndFishProcesses.z_euphotic; 
z_epipelagic = inputForPlots.carbonFluxAndFishProcesses.z_epipelagic; 

[figures.carbonFluxAndFishProcesses, ...
tilingLayouts.carbonFluxAndFishProcesses, ...
ax.carbonFluxAndFishProcesses] = plotInTiles(plotData,zSpan, ...
'Layout', {1,2}, ...
'FigureName','CarbonFluxAndFishProcesses', ...
'PlotOrder',1:2, ...
'LineSpecs',LineSpecs, ...
'Colors',Colors, ...
'Titles',Titles, ...
'XLabels',XLabels, ...
'YLabels',YLabels, ...
'DepthPlot',true, ...
'Legends',Legends);

% Adding helper lines:
xline(ax.carbonFluxAndFishProcesses(1), ...
    0,'Color',black,'LineStyle','--','LineWidth',1) % Zero flux.
yline(ax.carbonFluxAndFishProcesses(1), ...
    z_euphotic,'Color',black, ...
      'LineStyle','--','LineWidth',1) % Euphotic depth.
yline(ax.carbonFluxAndFishProcesses(1), ...
    z_epipelagic,'Color',black, ...
     'LineStyle','--','LineWidth',1) % End of epipelagic zone.

% Adding descriptive text to helper lines:
% Defining the x-coordinate position for the text boxes:
xPos = ax.carbonFluxAndFishProcesses(1).XTick(1) ...
    + 20*ax.carbonFluxAndFishProcesses(1).TickLength(1);

% Adding the text boxes: 
text(ax.carbonFluxAndFishProcesses(1), ...
    xPos,z_euphotic+10,'Euphotic depth','FontSize',9)
text(ax.carbonFluxAndFishProcesses(1), ...
    xPos,z_epipelagic+10,'200m','FontSize',9)

%% Set up for figure 4

% Defining options for the Vertical Carbon Flux plot:
Legends{1} = {'Relative gravity flux', 'Relative diffusive flux', ...
    'Relative active flux'};
Colors{1} = [red; 0.5*(yellow + green); blue]; 
LineSpecs{1} = {'-', '-', '-'};
Titles{1} = 'Vertical Relative Carbon Fluxes';
XLabels{1} = 'Relative downward flux (dimensionless)';
YLabels{1} = 'Depth (m)';

%% Figure 4: Vertical profile of relative carbon fluxes

plotData = inputForPlots.relativeCarbonFluxes.plotData; 
zSpan = inputForPlots.relativeCarbonFluxes.zspan; 
z_euphotic = inputForPlots.relativeCarbonFluxes.z_euphotic; 
z_epipelagic = inputForPlots.relativeCarbonFluxes.z_epipelagic; 

[figures.relativeCarbonFluxes, ...
tilingLayouts.relativeCarbonFluxes, ...
ax.relativeCarbonFluxes] = plotInTiles(plotData,zSpan, ...
'FigureName','RelativeCarbonFluxes', ...
'LineSpecs',LineSpecs, ...
'Colors',Colors, ...
'Titles',Titles, ...
'XLabels',XLabels, ...
'YLabels',YLabels, ...
'DepthPlot',true, ...
'Legends',Legends);

% Adding helper lines:
yline(ax.relativeCarbonFluxes(1), ...
    z_euphotic,'Color',black, ...
      'LineStyle','--','LineWidth',1) % Euphotic depth.
yline(ax.relativeCarbonFluxes(1), ...
    z_epipelagic,'Color',black, ...
     'LineStyle','--','LineWidth',1) % End of epipelagic zone.

% Adding descriptive text to helper lines:
% Defining the x-coordinate position for the text boxes:
xPos = ax.relativeCarbonFluxes(1).XTick(1) ...
    + 20*ax.relativeCarbonFluxes(1).TickLength(1);

% Adding the text boxes: 
text(ax.relativeCarbonFluxes(1), ...
    xPos,z_euphotic+10,'Euphotic depth','FontSize',9)
text(ax.relativeCarbonFluxes(1), ...
    xPos,z_epipelagic+10,'200m','FontSize',9)

%% Set up for figure 5 

%%% The irradiance:
XScale{1} = 'linear';
LineSpecs{1} = {'-','-'};  
Colors{1} = [blue; red];
Titles{1} = 'Light penetration (log-fraction)';
Legends{1} = {'Simulations','Observations'};

%%% The dissolved inorganic phosphorus DIP:
XScale{2} = 'linear';
LineSpecs{2} = {'-','--*'};  
Colors{2} = [blue; red];
Titles{2} = 'DIP (\mu M-P)';
Legends{2} = {'Simulations','Observations'}; 

%%% The silicate S:
XScale{3} = 'linear';
LineSpecs{3} = {'-','--*'}; 
Colors{3} = [blue; red];
Titles{3} = 'Silicate (\mu M-Si)';
Legends{3} = {'Simulations','Observations'}; 

%%% The dissolved organic carbon DOC:
XScale{4} = 'linear';
LineSpecs{4} = {'-','--*'}; 
Colors{4} = [blue; red];
Titles{4} = 'DOC (\mu M-C)';
Legends{4} = {'Simulations','Observations'}; 

%%% The dissolved oxygen:
XScale{5} = 'linear';
LineSpecs{5}  = {'-','-'};
Colors{5} = [blue; red];
Titles{5} = 'Dissolved oxygen (mL L^{-1})'; 
Legends{5} = {'Simulations','Observations'}; 

%%% The autotrophs
XScale{6} = 'linear';
LineSpecs{6} = {'-','-'};
Colors{6} = [blue; red];
Titles{6} = 'Autotrophs (rel. distr.)';
Legends{6} = {'Simulations','Observations'}; 

%%% The bacteria: 
XScale{7} = 'linear';
LineSpecs{7} = {'-','--*'};
Colors{7} = [blue; red];
Titles{7} = 'Bacteria (rel. distr.)';
Legends{7} = {'Simulations','Observations'}; 

%%% The mesozooplankton: 
XScale{8} = 'linear';
LineSpecs{8} = {'-','--*'};
Colors{8} = [blue; red];
Titles{8} = 'Mesozooplankton (rel. distr.)';
Legends{8} = {'Simulations','Scaled observations'}; 

%%% The gravitational POC flux:
XScale{9} = 'linear';
LineSpecs{9} = {'-','--','*'};
Colors{9} = [blue; blue; red];
Titles{9} = 'POC flux (mg C m^{-2}h^{-1})'; 
Legends{9} = {'Simulated partical gravity flux','Simulated particle flux', ...
                      'Observed'};

%%% The detritus:
XScale{10} = 'log';
LineSpecs{10} = {'-','*','-','*','-','*'};
Colors{10} = [red; red; blue; blue; green; green];
Titles{10} = 'Detritus composition (fraction)';
Legends{10} = {'Suspended','Observed','Slow sinking','Observed', ...
                        'Fast sinking','Observed'};

%% Figure 5: Observations versus simulations

x = inputForPlots.ObservationsAndSimulations.x; 
PlotData = inputForPlots.ObservationsAndSimulations.PlotData; 
gravityFlux = inputForPlots.ObservationsAndSimulations.gravityFlux; 
depths = inputForPlots.ObservationsAndSimulations.depths; 
standardDeviations ...
    = inputForPlots.ObservationsAndSimulations.standardDeviations;

[figures.ObservationsAndSimulations, ...
tilingLayouts.ObservationsAndSimulations, ...
ax.ObservationsAndSimulations] = plotInTiles(PlotData,x, ...
'Layout', {2,5}, ...
'FigureName','ObservationsAndSimulations', ...
'XScale',XScale, ...
'LineSpecs',LineSpecs, ...
'Colors',Colors, ...
'Titles',Titles, ...
'DepthPlot',true, ...
'Legends',Legends);

% Adding error bar to the POC plot: 
hold(ax.ObservationsAndSimulations(9),"on")
errorbar(ax.ObservationsAndSimulations(9),gravityFlux,depths, ...
    gravityFlux-standardDeviations,gravityFlux+standardDeviations, ...
    'horizontal','*','Color',red)
ax.ObservationsAndSimulations(9).Legend.Location = 'best';
hold(ax.ObservationsAndSimulations(9),"off")
end
