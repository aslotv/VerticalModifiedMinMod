function [tiledPlot,tilingLayout,ax] = plotInTiles(x,y,options)
%PLOTINTILES is a function for making a figure with tiled plots. There are
%two required inputs, x and y, which must be given as 1 by N cell arrays,
%where N is the number of plots in the tiled layout. The x and y data in
%each cell, must be given as 1 by J cells, where J is the number of
%line-plots in the tiled plot, and where the contents of each of the J
%cells must be vectors. This functions stores and returns the figure
%handle, the tiled layout handle, as well as the axis handle.
%
%   Inputs:
%
%   Required:
%   x - a 1 by N cell array, where N is the number of plots in
%   the tiled layout, of inputs for the x-axis in the plot function stored
%   in a cell.
%
%   y - a 1 by N cell array of inputs for the y-axis in the plot
%   function, stored in a cell.
%
%   Optional:
%   options - A set of name-value argument pairs, which can be
%                  given as either name=value, or 'name',value.
%
%   Name-Value pairs:
%   FigureName (No name as default) - Must be given as a text scalar.
%
%   MainTitle (No title as default) - Must be given as a text scalar.
%
%   MainXLabel (No label as default) - Must be given as a text scalar.
%
%   MainYLabel (No label as default) - Must be given as a text scalar.
%
%   Layout ('flow' as default) - Must be either 'flow' or a two dimensional
%   row vector [m n], where m is the number of rows in the layout grid,
%   and n is the number of columns (m*n must equal N).
%
%   PlotOrder (No re-ordereing as default) - A row-vector of length N, with
%   integers from 1 to N.
%
%   XScale ('linear' as default) - A  1 by N cell of text scalars; must be
%   either 'linear' or 'log'.
%
%   YScale ('linear' as default) - A  1 by N cell of text scalars; must be
%   either 'linear' or 'log'.
%
%   Titles (No titles as default) - A  1 by N cell of text scalars.
%
%   XLabels (No labels as default) - A  1 by N cell of text scalars.
%
%   YLabels (No labels as default) - A  1 by N cell of text scalars.
%
%   Legends (No labels as default) - A 1 by N cell, where each cell
%   contains a cell of text scalars.
%
%   LineSpecs (Continuous line as default) - A 1 by N cell, where each cell
%   contains a cell of line specs (e.g., '-', '-.',  '--*', or '*') for
%   each line in plot number i, i=1,2,...,N.
%
%   LineWidth (LineWidth = 1  as default) - Must be a scalar.
%
%   Colors (MATLAB(R) default plot colors as default) - Must be given as a
%   1 by N cell, where each cell contains an m by 3 matrix of RGB triplets,
%   where m is the number of colors.
%
%   DepthPlot (DepthPlot = false as default) - Must be given as a single
%   logical scalar (true or false), or as a 1 by N cell of logicals, where
%   each cell corresponds to the inputs x and y.
%
%   Output:
%
%   tiledPlot - The figure handle (same as f = figure).
%
%   tilingLayout - The tiled layout handle (same as t = tiledlayout).
%   
%   ax - The axis handels for each tiled plot (same as ax(n) = nexttile).
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

%% Defining function arguments
arguments
    %%% Required input:
    x (1,:) {mustBeA(x,'cell')}
    y (1,:) {mustBeA(y,'cell')}

    %%% Optional name-value argument input:

    % Name of figure, and main title and axis labels:
    options.FigureName {mustBeTextScalar}
    options.MainTitle {mustBeTextScalar}
    options.MainXLabel {mustBeTextScalar}
    options.MainYLabel {mustBeTextScalar}

    % The order and layout of the tiled plots:
    options.Layout {mustBeA(options.Layout,'cell')}
    options.PlotOrder (1,:) {mustBeInteger}

    % The type of plot:
    options.XScale (1,:) {mustBeText}
    options.YScale (1,:) {mustBeText}

    % Titles and axis labels, and legends of tiled plots:
    options.Titles (1,:) {mustBeText,mustBeA(options.Titles,'cell')}
    options.XLabels (1,:) {mustBeText,mustBeA(options.XLabels,'cell')}
    options.YLabels (1,:) {mustBeText,mustBeA(options.YLabels,'cell')}
    options.Legends {mustBeA(options.Legends,'cell')}

    % Appearance of tiled plots:
    options.LineSpecs (1,:) {mustBeA(options.LineSpecs,'cell')}
    options.LineWidth (1,1) {mustBePositive}
    options.Colors {mustBeA(options.Colors,'cell')}
    options.DepthPlot (1,:) {mustBeA(options.DepthPlot,["logical","cell"])}

end

% Number of tiled plots in tiled layout:
NumPlots = length(y);

% Assigning default values to missing fields in options:
options = setMissingFields(options);

% Defining options.DepthPlot as cells, each containing value of DepthPlot,
% if DepthPlot is given as singal logical value:
if ~iscell(options.DepthPlot)
    DepthPlot = options.DepthPlot;
    options.DepthPlot = cell(1,NumPlots);
    [options.DepthPlot{1:NumPlots}] = deal(DepthPlot);
end


%% Creating the tiled plot

% Creating a new figure, with FigureName as name and setting the window
% state to maximized, and storing it in the handle tiledPlot:
tiledPlot = figure('Name',options.FigureName,'WindowState',"maximized");

pause(0.5) % Allowing figure properties time to settle.

% Initiates a tiled layout and stores the handle:
tilingLayout = tiledlayout(tiledPlot,options.Layout{:}, ...
    "TileSpacing","tight", ...
    "Padding","compact");

% Creating the tiled plots:
for n = 1:length(options.PlotOrder)

    % The index k is the index at the j-th position in the vector
    % PlotOrder:
    k = options.PlotOrder(n);

    if iscell(x{k}) && iscell(y{k}) && ~iscell(options.LineSpecs{k})
        [options.LineSpecs{1:length(x{k})}] = deal(options.LineSpecs(k));
    end

    % Putting the input for the plot function in a comma separated list:
    plotCell = reshape([x{k}; y{k}; options.LineSpecs{k}],1,[]);

    % Initiating a new tile in tilingLayout and storing the axis handle:
    ax(n) = nexttile; %#ok<*AGROW>

    % Making the j-th plot:
    plot(ax(n),plotCell{:},'LineWidth',options.LineWidth)

    % Defining axis-scales:
    set(gca,'XScale',options.XScale{k},'YScale',options.YScale{k})

    % Adding a title, and axis labels for the tiled plot:
    title(options.Titles{k},'FontWeight','bold')
    xlabel(options.XLabels{k},'FontWeight','bold')
    ylabel(options.YLabels{k},'FontWeight','bold')

    grid on

    % Adjusting the axes and limits:
    if options.DepthPlot{k}

        % Flipping the y-axis:
        axis ij

        % Finding the limits of the y axis:
        yMax = max(cellfun(@max,y{k}));
        yMin = min(cellfun(@min,y{k}));

        % Estimating the size order of the upper limit:
        szOrder = round(log10(yMax)) - 1;

        % Defining rounded limits:
        yMax = ceil(yMax/10^szOrder)*10^szOrder;
        yMin = floor(yMin/10^szOrder)*10^szOrder;

        if ~isequal(yMin,yMax)

            % Finding the number of digits in the range of the y limits:
            tickStep = 10^(round(log10(abs(yMax-yMin)))-1);

            % Defining the y limits and y ticks:
            ylim([yMin yMax])
            yticks(fliplr(yMax:-tickStep:yMin))

            % Defining ticklabels on every other tick:
            N_ticks = length(fliplr(yMax:-tickStep:yMin)); % # of ticks.
            ticklabels = strings(1,N_ticks); % Empty string array.
            ticklabels(fliplr(N_ticks:-2:1)) ...
                = string(fliplr(yMax:-2*tickStep:yMin)); % Adda tick labels
            ticklabels = num2cell(ticklabels); % Storing in cell array.
            yticklabels(ticklabels) % Assigning ticklabels to y-axis.

        end

        % The range of the x axis limits:
        xMax = max(cellfun(@max,x{k})); % The maximum x axis limit.
        xMin = min(cellfun(@min,x{k})); % The minimum x axis limit.
        range = xMax-xMin; % The range of the x axis limits.

        % Specifying the mode of the x limits as padded if the range in the
        % values of the x axis is larger than smallest positive normalized
        % floating point number of class double, if not the default mode
        % stands:
        if range > realmin
            xlim padded
        end

    else
        xlim padded

        % The range of the x axis limits:
        yMax = max(cellfun(@max,y{k})); % The maximum x axis limit.
        yMin = min(cellfun(@min,y{k})); % The minimum x axis limit.
        range = yMax-yMin; % The range of the x axis limits.

        % Specifying the mode of the x limits as padded if the range in
        % the values of the x axis is larger than smallest positive
        % normalized floating point number of class double, if not the
        % default mode stands: 
        if range > realmin
            ylim padded
        end
    end

    % Sets legends if legends are specified:
    if any(cellfun(@(x) x~="",options.Legends{k}))
        legend(options.Legends{k},'Location','best','AutoUpdate','off')
    end

    colororder(ax(n),options.Colors{k})
end

% Adding main title and axis labels:
title(tilingLayout,options.MainTitle,'FontWeight','bold')
xlabel(tilingLayout,options.MainXLabel,'FontWeight','bold')
ylabel(tilingLayout,options.MainYLabel,'FontWeight','bold')

%% Nested functions
    function options = setMissingFields(options)
        % This function finds the options that wasn't specified as name-value
        % argument pairs, and sets them to the default value.

        % The names of all the fields in the options structure:
        fieldNames = split("FigureName MainTitle " + ...
            "MainXLabel MainYLabel Layout PlotOrder  XScale YScale " + ...
            "Titles XLabels YLabels Legends LineSpecs " + ...
            "LineWidth Colors DepthPlot");

        % No figure name, main title or axis labels as default:
        [default.FigureName, default.MainTitle, ...
            default.MainXLabel, default.MainYLabel] = deal("");

        % The default layout setting for the tiling is 'flow':
        default.Layout = {'flow'};

        % No re-ordering as default of tiled plots as default:
        default.PlotOrder = 1:NumPlots;

        % The axes-scales are set to linear as default:
        [default.XScale{1:NumPlots}] = deal('linear');
        [default.YScale{1:NumPlots}] = deal('linear');

        % No titles or axis labels as default:
        [default.Titles{1:NumPlots}] = deal("");
        [default.XLabels{1:NumPlots}] = deal("");
        [default.YLabels{1:NumPlots}] = deal("");
        [default.Legends{1:NumPlots}] = deal("");

        % Continuous line as default:
        default.Linespecs = cell(1,NumPlots);
        for i = 1:NumPlots
            [default.LineSpecs{i}{1:length(x{i})}] = deal('-');
        end

        % Using double of MATLAB(R) default LineWidth (0.5) as default:
        default.LineWidth = 1;

        % Using MATLAB(R) default colororder as default:
        [default.Colors{1:NumPlots}] = deal('default');

        % No DepthPlot as default:
        [default.DepthPlot{1:NumPlots}] = deal(false);

        % Assigns default values for the missing fields in the options
        % structure:
        for i = 1:length(fieldNames)
            if ~isfield(options,fieldNames(i))
                options.(fieldNames(i)) = default.(fieldNames(i));
            end
        end

    end

pause(0.5) % Allowing figure properties time to settle.
end
