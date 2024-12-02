% ===========================================================
% File: exportFigures.m
% Description: This script exports the figures based on simulations for the
%              mean squared error, bias, variance, weights, and
%              unrestricted units for unit averaging.
%
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Plot creation is split into three sections:
%   I. Plots reported in the main text: relative mean squared error, 
%      average weight assigned to the target unit, side-by-side comparison 
%      of bias and variance. 
%  II. Line plots reported in the online appendix: mean squared error,
%      full results for bias and variance. This section includes
%      comparisons for various choices of tuning parameters in large-N
%      approaches
% III. Grid plots reported in the online appendix: unit weight as a
%      function of own parameter value and target parameter value;
%      probability of being unrestricted as a function of own and target
%      value.
%  IV. An animated plot for discussing results informally. Combines
%      information on the MSE, bias, and variance.
%
% Created figures are saved in the 'results/figures/' folder.
% ===========================================================

%% Initialization of plotting

% Load current file
fileSaveName = makeOutputFileName('MSE', simulationSetting, ...
                                  numReplicationsMSE, valuesN, valuesT);
load(fileSaveName)

% Load bundles of averaging approaches used in plotting
setApproachesToPlot

% Set parameters for image creation
setPlottingParameters

% Load line and grip plot descriptions
chooseLinePlots
chooseGridPlots
   
%% Main Text Plots: MSE and average first weight

% Set line plots to use lines for the paper
linePlots = changeLinesToPlot(linePlots, approachesToPlotMSEPaper);

% Main text: report the results only for lambda_1
parID = 1;

% --- MSE and average first weight plots ---

% IDs of plots to create in this plot
textPlotIDs = {1, 5};

% Values of N and T to use
valuesNPlotArray{1} = logical([1, 1, 1]);
valuesTPlotArray{1} = logical([1, 1, 0, 1]);

valuesNPlotArray{2} = logical([1, 1, 0]);
valuesTPlotArray{2} = logical([0, 1, 0, 1]);

% Plot saving sizes
textPlotSavingSize{1} = [0 0 12 11]*0.9;
textPlotSavingSize{2} = [0 0 22 10]*0.5;

% Upper margin
textPlotUpperMargin{1} = 0.12;
textPlotUpperMargin{2} = 0.16;

% Loop through the plots to be created 
for textPlotID = 1:length(textPlotIDs)
    % Extract plot ID and value of (N, T) to use
    plotID = textPlotIDs{textPlotID};
    valueNPlot = valuesNPlotArray{textPlotID};
    valueTPlot = valuesTPlotArray{textPlotID};
    numTPlot = sum(valueTPlot);
    numNPlot = sum(valueNPlot);
    
    % Extract description of the plot
    plotDescrStruct = linePlots{plotID};
    
    % Extract the short names of the approaches for this plot
    approachesToPlot = plotDescrStruct.approachesToPlot;
    
    % Obtain parameter name
    paramName = paramArray{parID}.saveName;
    
    % Create figure
    if plotQuietly ~= 1
        figure('Renderer', 'painters', ...
            'Position', [50 50 1900 1200]); %#ok<*FGREN>
    else
        figure('visible','off',...
            'Renderer', 'painters',...
            'Position', [50 50 plotWLines plotHLines]);
    end
    
    % Create tight axes using above parameters
    [ha, ~] = tight_subplot(numTPlot,numNPlot, ...
        [.052 .023],[.07 textPlotUpperMargin{textPlotID}],[.05 .01]);
    
    % Initialize counters for values of (N, T)
    tPlotID = 0;
    nPlotID = 0;

    % Loop over N and T included in this plot
    for tID = 1:numT % T indexes rows
        
        % Skip iteration is current value of T is not required
        if ~valueTPlot(tID)
            continue
        end
        tPlotID = tPlotID + 1;
        
        for nID = 1:numN % N indexes columns
            
            % Skip iteration is current value of N is not required
            if ~valueNPlot(nID)
                continue
            end
            
            % Compute number of plot and select corresponding axes
            nPlotID = mod(nPlotID, numNPlot)+1;
            plotNum = (tPlotID-1)*numNPlot+nPlotID;
            axes(ha(plotNum)) %#ok<*LAXES>
            
            % Extract data
            currentParTable = plotDescrStruct.data{nID, tID}.(paramName);
            
            % Extract approaches present in this table
            approachesPresent = currentParTable.Properties.VariableNames;
            
            % Extract number of approaches
            numApproaches = size(currentParTable, 2);
            
            % Set hold on
            hold(ha(plotNum), 'on');
            
            % Loop through the columns currentParTable (approaches present)
            for approachID = plotDescrStruct.firstColumnPlot:numApproaches
                
                % Get approach short name
                approachShortName = string(approachesPresent{approachID});
                
                % Check if approach is present in the list of
                % approaches to plot
                approachPresent = ismember(approachShortName,  ...
                    approachesToPlot);
                % Skip if this approach was not available or not present
                if ~approachPresent
                    continue
                end
                
                % Find coordinate of current column/approach in the
                % description of the approaches
                approachCoord = findApproach(allMethodsArray, ...
                    approachShortName);
                % Extract approach struct
                approach = allMethodsArray{approachCoord};
                % Obtain name
                approachName = approachPlotName(approach, ...
                    approachPaperRenames, ...
                    approachShortName);
                
                % Obtain data to plot according to the descriptions of the
                % plot
                currentLineData = ...
                    plotDescrStruct.dataTransform(...
                    currentParTable, approachID);
                
                % Smooth out the line (only affects the odd edge in the
                % mean group and AIC curves for the MSE)
                dataInterp = ...
                    interp1(theta1Range, currentLineData, ...
                    thetaGridMSE, 'spline');
                
                % Plot the smoothed out line
                hL1 = plot(thetaGridMSE, dataInterp, ...
                    'LineStyle', approach.lineStyle, ...
                    'LineWidth', plotLineThickness, ...
                    'Color', approach.colorBW);
                hL1.HandleVisibility = 'off';
                
                % Add the markers
                hL2 =  plot(theta1Range, currentLineData, ...
                    'LineStyle', 'none', ...
                    'Marker', approach.marker, ...
                    'MarkerSize', approach.markerSize, ...
                    'Color', [0,0,0],... %approach.colorBW, ...
                    'DisplayName', approachName);
                hL2.HandleVisibility = 'off';
                
                % Dummy line with both line and markers for the legend
                hLDummy = plot(2, 2, ...
                    'LineStyle', approach.lineStyle, ...
                    'LineWidth', plotLineThickness, ...
                    'Marker', approach.marker, ...
                    'Color', approach.colorBW, ...
                    'DisplayName', approachName);
            end
            
            % Set limits for axes
            xlim([min(theta1Range), max(theta1Range)])
            yLims = plotDescrStruct.yLims(currentParTable, paramName);
            ylim(yLims)
            
            % Add vertical line if plot is relative
            if plotDescrStruct.relative
                hV = yline(1);
                hV.HandleVisibility='off';
            end
            
            % Add a title and left-justify it
            ttl = title("N="+num2str(valuesN(nID))+...
                ", T="+num2str(valuesT(tID)));
            ttl.Units = 'Normalize';
            ttl.Position(1) = 0;
            ttl.HorizontalAlignment = 'left';
            
            % Add x axis labels only on the bottom
            if tID == numT
                xlabel('$\lambda_1$')
                ha(plotNum).XTickLabel = ha(plotNum).XTick;
            end

            % Add y axis labels only on the left
            if nID == 1
                yTicks = determineTicks(yLims, 4, 10);
                ha(plotNum).YTick = yTicks;
                ha(plotNum).YTickLabel = yTicks;
            end

            % Enforce same y-axis ticks on all plots
            ha(plotNum).YTick = yTicks;
        end
    end
    
    % Add the suptitle to the plot
    sgtitle(plotDescrStruct.suptitle(paramArray{parID}.plotDescr))
    
    % Add legend
    legend('Position', [0.83, 0.88, 0.1, 0.06]);
    
    % Create a name to the export the figure
    figureSavingName =  makeFigureFileName(...
        'BW', simulationSetting, ...
        plotDescrStruct.plotSavingName, paramName, ...
        valuesN, valuesT);
    
    % Save the figure in PNG 
    print(gcf, figureSavingName, '-dpng', '-r300' );

    % Save the figure in PDF
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto',...
        'PaperUnits','Inches',...
            'PaperSize',[pos(3), pos(4)])
    print(gcf, figureSavingName, '-dpdf' );
end

% --- Bias-variance plot ---
 
% Set value of (N, T) to use
nID = 2; % N=150
tID = 2; % T=60

% Create figure
if plotQuietly ~= 1
    figure('Renderer', 'painters', ...
        'Position', [50 50 1200 800]);
else
    figure('visible','off',...
        'Renderer', 'painters',...
        'Position', [50 50 plotWLines plotHLines]);
end

% Create tight axes
[ha, ~] = tight_subplot(1, 2, ...
    [.04 .04],[.14 0.08],[.05 .01]);

% plotID 3 and 4 correspond to bias and variance
for plotID = 3:4
    % Axes ID to use
    plotNum = plotID-2;
    axes(ha(plotNum))  
    
    % Extract description of the plot
    plotDescrStruct = linePlots{plotID};
    
    % Extract the approaches for this plot
    approachesToPlot = plotDescrStruct.approachesToPlot;
    
    % Extract data
    currentParTable = plotDescrStruct.data{nID, tID}.(paramName);
    
    % Extract approaches present in this table
    approachesPresent = currentParTable.Properties.VariableNames;
    
    % Extract number of approaches
    numApproaches = size(currentParTable, 2);
    
    % Plot
    hold(ha(plotNum), 'on');
    % Loop through the columns of the currentParTable
    for approachID = plotDescrStruct.firstColumnPlot:numApproaches
        
        % Get approach short name
        approachShortName = string(approachesPresent{approachID});
        
        % Check if approach is present in the list of
        % approaches to plot
        approachPresent = ismember(approachShortName ,approachesToPlot);
        % Skip if this approach was not available or not required
        if ~approachPresent
            continue
        end
        
        % Find coordinate of current column/approach in the
        % description of the approaches
        approachCoord = findApproach(allMethodsArray, ...
            approachShortName);
        % Extract approach struct
        approach = allMethodsArray{approachCoord};
        % Obtain name
        approachName = approachPlotName(approach, ...
            approachPaperRenames, ...
            approachShortName);
        
        % Obtain data to plot according to the descriptions of the
        % plot
        currentLineData = ...
            plotDescrStruct.dataTransform(...
            currentParTable, approachID);
        
        % Interpolate data between markers
        dataInterp = ...
            interp1(theta1Range, currentLineData, ...
            thetaGridMSE, 'spline');
        
        % Plot the smoothed data line
        hL1 = plot(thetaGridMSE, dataInterp, ...
            'LineStyle', approach.lineStyle, ...
            'LineWidth', plotLineThickness + ...
            1*(approachName=="Individual"), ...
            'Color', approach.colorBW);
        hL1.HandleVisibility = 'off';
        
        % Add the markers
        hL2 =  plot(theta1Range, currentLineData, ...
            'LineStyle', 'none', ...
            'Marker', approach.marker, ...
            'MarkerSize', approach.markerSize, ...
            'Color', [0,0,0],... %approach.colorBW, ...
            'DisplayName', approachName);
        hL2.HandleVisibility = 'off';
        
        % Dummy line with both markers and line, for legend
        hLDummy = plot(2, 2, ...
            'LineStyle', approach.lineStyle, ...
            'LineWidth', plotLineThickness, ...
            'Marker', approach.marker, ...
            'Color', approach.colorBW, ...
            'DisplayName', approachName);
        
        
    end
    % Set limits for axes
    xlim([min(theta1Range), max(theta1Range)])
    yLims = plotDescrStruct.yLims(currentParTable, paramName);
    ylim(yLims)
    
    % Add vertical line if appropriate
    if plotDescrStruct.relative
        hV = yline(1);
        hV.HandleVisibility='off';
    end
    
    % Add special titles
    if plotID == 3
        ttl = title("Bias");
    elseif plotID == 4
        ttl = title("Variance relative to individual estimator");
    end
    % Left-justify the title
    ttl.Units = 'Normalize';
    ttl.Position(1) = 0;
    ttl.HorizontalAlignment = 'left';
    
    % Add x label and x ticks
    xlabel('$\lambda_1$')
    ha(plotNum).XTickLabel = ha(plotNum).XTick;
    
    % Determine y ticks separately for each plot
    yTicks = determineTicks(yLims, 3, 10);
    ha(plotNum).YTick = yTicks;
    ha(plotNum).YTickLabel = yTicks;
    
    % Add legend on the bias pane (to include ind estimator)
    if plotNum ==1
        legend('Position', [0.34, 0.76, 0.1, 0.06]);
    end
end

% Create a name to the export the figure
figureSavingName =  makeFigureFileName(...
    'BWbias_variance', simulationSetting, ...
    plotDescrStruct.plotSavingName, paramName, ...
    valuesN, valuesT);

% Save the figure in PNG
print(gcf, figureSavingName, '-dpng', '-r300' );

% Save the figure in PDF
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto',...
    'PaperUnits','Inches',...
    'PaperSize',[pos(3), pos(4)])
print(gcf, figureSavingName, '-dpdf' );

%% Online Appendix: Line Plots

% Loop through sets of averaging approaches
for lineSetID = 1:length(lineSets)
    
    % Set line plots to use
    linePlots = changeLinesToPlot(linePlots, lineSets{lineSetID});
    
    % Loop through plot descriptions
    for plotID = 1:5

        % Extract description of the plot
        plotDescrStruct = linePlots{plotID};
        
        % Extract the approaches for this plot
        approachesToPlot = plotDescrStruct.approachesToPlot;
        
        % Loop over parameter   s
        for parID = 1:numParams

            % For main set of lines, generate MSE and weights for all, and
            % bias and variance only for lambda
            if parID > 1 && ismember(plotID, [2, 3, 4]) 
                continue
            end
 
            % Parameter name
            paramName = paramArray{parID}.saveName;
            
            % Create figure
            if plotQuietly ~= 1
                figure('Renderer', 'painters', ...
                    'Position', [50 50 plotWLines plotHLines]);
            else
                figure('visible','off',...
                    'Renderer', 'painters',...
                    'Position', [50 50 plotWLines plotHLines]);
            end
            
            [ha, ~] = tight_subplot(numT,numN,[.04 .023],[.07 .01],[.05 .01]);
            
            % Loop over N and T
            for tID = 1:numT
                for nID = 1:numN
                    % Create subplot
                    plotNum = (tID-1)*numN+nID;
                    axes(ha(plotNum)) %, 'Parent', p);
                    
                    % Extract data
                    currentParTable = plotDescrStruct.data{nID, tID}.(paramName);
                    
                    % Extract approaches present in this table
                    approachesPresent = ...
                        currentParTable.Properties.VariableNames;
                    
                    % Extract number of approaches
                    numApproaches = size(currentParTable, 2);
                    
                    % Plot
                    hold(ha(plotNum), 'on');
                    for approachID = plotDescrStruct.firstColumnPlot:numApproaches
                        
                        % Get approach short name
                        approachShortName = ...
                            string(approachesPresent{approachID});
                        
                        % Check if approach is present in the list of
                        % approaches to plot
                        approachPresent = ...
                            ismember(approachShortName ,approachesToPlot);
                        % Skip if this approach was not available
                        if ~approachPresent
                            continue
                        end
                        
                        % Find approach and extract information about it
                        approachCoord = ...
                            findApproach(allMethodsArray, approachShortName);
                        approach = allMethodsArray{approachCoord};
                        approachName = approachPlotName(approach, ...
                            nameSets{lineSetID}, ...
                            approachShortName);
                        
                        % Compute data for the current approach
                        currentLineData = ...
                            plotDescrStruct.dataTransform(...
                            currentParTable, approachID);
                        
                        % Smooth out the line (only affects the odd edge in  
                        % the mean group and AIC curves for the MSE)
                        dataInterp = ...
                            interp1(theta1Range, currentLineData, ...
                            thetaGridMSE, 'spline');

                        % Plot the smoothed out line
                        hL1 = plot(thetaGridMSE, dataInterp, ...
                            'LineStyle', approach.lineStyle, ...
                            'LineWidth', plotLineThickness, ...
                            'Color', approach.color);
                        hL1.HandleVisibility = 'off';
                        
                        % Add the markers
                        hL2 =  plot(theta1Range, currentLineData, ...
                            'LineStyle', 'none', ...
                            'Marker', approach.marker, ...
                            'MarkerSize', approach.markerSize, ...
                            'Color', approach.color, ...
                            'DisplayName', approachName);
                        hL2.HandleVisibility = 'off';
                        
                        % Dummy line for legend
                        hLDummy = plot(2, 2, ...
                            'LineStyle', approach.lineStyle, ...
                            'LineWidth', plotLineThickness, ...
                            'Marker', approach.marker, ...
                            'Color', approach.color, ...
                            'DisplayName', approachName);
                    end 

                    % Set limits for axes
                    xlim([min(theta1Range), max(theta1Range)])
                    if nID == 1
                        yLims = plotDescrStruct.yLims(...
                            currentParTable, paramName);
                    end
                    ylim(yLims)
                    
                    % Add vertical line if appropriate
                    if plotDescrStruct.relative
                        hV = yline(1);
                        hV.HandleVisibility='off';
                    end 
                    
                    % Add a title and left-justify it
                    ttl = ...
                        title("N="+num2str(valuesN(nID))+ ...
                            ", T="+num2str(valuesT(tID)));
                    ttl.Units = 'Normalize';
                    ttl.Position(1) = 0;  
                    ttl.HorizontalAlignment = 'left';
                    
                    % x axis labels only on the bottom
                    if tID == numT
                        xlabel('$\lambda_1$')
                        ha(plotNum).XTickLabel = ha(plotNum).XTick;
                    end
                    % y axis labels only on the left
                    if nID == 1 
                        yTicks = determineTicks(yLims, 3, 10);
                        ha(plotNum).YTick = yTicks;
                        ha(plotNum).YTickLabel = yTicks;
                    end
                    % Enforce same ticks in all subplots
                    ha(plotNum).YTick = yTicks;
                end
            end
            
            % Include suptitle
            sgtitle(plotDescrStruct.suptitle(paramArray{parID}.plotDescr))
            
            % Add legend
            legend('Position', [0.845, 0.9, 0.08, 0.066]);

            % Create a name to the export the figure
            figureSavingName =  makeFigureFileName(...
                plotName{lineSetID}, simulationSetting, ...
                plotDescrStruct.plotSavingName, paramName, ...
                valuesN, valuesT);

            % Save the figure in PNG
            print(gcf, figureSavingName, '-dpng', '-r300' );

            % Save the figure in PDF
            set(gcf,'Units','Inches');
            pos = get(gcf,'Position');
            set(gcf,'PaperPositionMode','Auto',...
                'PaperUnits','Inches',...
                'PaperSize',[pos(3), pos(4)])
            print(gcf, figureSavingName, '-dpdf' );

        end
    end
end

%% Online Appendix: Grid Plots of Weights
 
% Grid results plotted only for the AR parameter
parID = 1;

% Loop through types
for plotID = 1:length(gridPlots)

    % Loop through averaging approaches
    for approachID = 1:length(gridPlots{plotID}.approachesToPlot)

        % Extract current averaging approach
        approachWeightPlot = gridPlots{plotID}.approachesToPlot(approachID);

        % Find coordinate of the current approache
        for candApproachID = 1:length(allMethodsArray)
            if allMethodsArray{candApproachID}.shortName == approachWeightPlot
                break
            end
        end

        % Extract the approach from the methods array and prepare its name
        approach = allMethodsArray{candApproachID};
        approachName = approachPlotName(approach, ...
            approachGridRenames, ...
            approachWeightPlot);

        % Extract parameter name
        paramName = paramArray{parID}.saveName;

        % Create figure
        if plotQuietly ~= 1
            p = figure('Renderer', 'painters', ...
                'Position', [50 50 plotWLines plotHLines]);
        else
            p = figure('visible','off',...
                'Renderer', 'painters',...
                'Position', [50 50 plotWLines plotHLines]);
        end
        % Create meshgrid
        [X, Y] = meshgrid(theta1Range, thetaGrid);

        % Create plotting grid
        [ha, ~] = tight_subplot(numT,numN,[.04 .023],[.07 .12],[.05 .01]);

        % Loop over N and T
        for tID = 1:numT
            for nID = 1:numN
                % Extract weights
                dataNW = ...
                    gridPlots{plotID}.data{nID, tID}.(paramName).(approachWeightPlot);

                % Create subplot
                plotNum = (tID-1)*numN+nID;
                axes(ha(plotNum))

                % Draw the results as a pseudocolor plot    
                pcolor(X, Y, dataNW) 

                % White for smallest values, blue for intermediate, orange 
                % for largest values
                colormap(...
                    orangeBlueMaps(...
                    gridPlots{plotID}.cmapN1,...
                    gridPlots{plotID}.cmapN2,...
                    gridPlots{plotID}.cmapFlip))

                % Set color limits
                clim(gridPlots{plotID}.caxisLims)  

                % Adjust the aspect ratio of the plots
                tmpAspect=daspect(); 
                daspect(tmpAspect([1 2 2]));

                % Change scaling to logarithmic if required
                if gridPlots{plotID}.logScale
                    set(gca,'ColorScale','log')
                end

                % Add a colobar on the right
                if nID == numN
                    colorbar
                end

                % Add a title and left-justify it
                ttl = title("N="+num2str(valuesN(nID))+...
                    ", T="+num2str(valuesT(tID)));
                ttl.Units = 'Normalize';
                ttl.Position(1) = 0;
                ttl.HorizontalAlignment = 'left';

                % Add axis labels on the outside plots
                if nID == 1
                    ylabel('$\lambda_1$')
                end
                if tID == numT
                    xlabel('$\lambda_{alt}$')
                end
            end
        end
        % Add suptitle
        sgtitle({gridPlots{plotID}.suptitle, ...
            approachName + ", " + paramArray{parID}.plotDescr})

        % Create a name to the export the figure
        figureSavingName =  makeFigureFileName(...
            approachWeightPlot, simulationSetting, ...
            gridPlots{plotID}.plotSavingName, paramName, ...
            valuesN, valuesT);

        % Save the figure in PNG
        print(gcf, figureSavingName, '-dpng', '-r300' );

        % Save the figure in PDF
        set(gcf,'Units','Inches');
        pos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto',...
            'PaperUnits','Inches',...
            'PaperSize',[pos(3), pos(4)])
        print(gcf, figureSavingName, '-dpdf' );
    end
end

%% Animated plot for discussing the simulation results

% Report animation for the AR(1) parameter
parID = 1;

% Set line plots to use lines for the animated plot
linePlots = changeLinesToPlot(linePlots, approachesToPlotAnimated);

% Report plots for MSE, bias, variance
plotCoords = [1, 3, 4];

% Sample sizes reported: (N, T) = (150, 60)
nID = 2;
tID = 2;

% Create figure
if plotQuietly ~= 1
    figure('Renderer', 'painters', ...
        'Position', [50 50 1200 420]); %#ok<*FGREN>
else
    figure('visible','off',...
        'Renderer', 'painters',...
        'Position', [50 50 plotWLines plotHLines]);
end

% Set figure background to white 
set(gcf, 'Color', 'w');  

% Create tight axes using above parameters
[ha, ~] = tight_subplot(1,3, ...
    [.052 .053],[.09 0.12],[.05 .01]);

% Create arrays for plot data. These will be used in the dynamic component
hLines = cell(3, 1);
hMarkers = cell(3, 1);
dataInterp = cell(3, 1);
currentLineData = cell(3, 1);

% Plot the static background
for plotID = 1 : 3
    % Extract the kind of plot desired
    plotKind = plotCoords(plotID);

    % Extract description of the plot
    plotDescrStruct = linePlots{plotKind};

    % Extract the short names of the approaches for this plot
    approachesToPlot = plotDescrStruct.approachesToPlot;

    % Obtain parameter name
    paramName = paramArray{parID}.saveName;
    % Extract data
    currentParTable = plotDescrStruct.data{nID, tID}.(paramName);

    % Extract approaches present in this table
    approachesPresent = currentParTable.Properties.VariableNames;

    % Extract number of approaches
    numApproaches = size(currentParTable, 2);

    % Switch axes
    axes(ha(plotID))

    % Set hold on
    hold(ha(plotID), 'on');

    % Add a box for nicer visuals
    box on

    % Statically add the individual and the mean group estimators
    % Also add an empty line for the fixed-N estimator so that it appears
    % on the legend.
    for approachID = ...
            plotDescrStruct.firstColumnPlot : length(approachesToPlot) 

        % Extract current approach
        approachShortName = approachesToPlot(approachID);

        % Find coordinate of current column/approach in the
        % description of the approaches
        approachCoord = findApproach(allMethodsArray, ...
            approachShortName);
        approachColumnID = ...
            find(string(currentParTable.Properties.VariableNames) == ...
            approachShortName);

        % Extract approach struct
        approach = allMethodsArray{approachCoord};

        % Obtain name
        approachName = approachPlotName(approach, ...
            approachPaperRenames, ...
            approachShortName);

        % Obtain data to plot according to the descriptions of the
        % plot
        currentLineData{plotID} = ...
            plotDescrStruct.dataTransform(...
            currentParTable, approachColumnID);

        % Smooth out the mean group line
        dataInterp{plotID} = ...
            interp1(theta1Range, currentLineData{plotID}, ...
            thetaGridMSE, 'spline');
        
        % Fixed-N: only dummy data at this stage
        if approachID == 3 
            dataLine = nan(length(thetaGridMSE), 1);
            dataMarker = nan(length(theta1Range), 1);
            plotLineThickness = 1.6;
        else
            % Other approaches: full data
            dataLine = dataInterp{plotID};
            dataMarker = currentLineData{plotID};
            plotLineThickness = 0.8;
        end    
        
        % Plot the smoothed out line
        hLines{plotID} = plot(thetaGridMSE, dataLine, ...
            'LineStyle', approach.lineStyle, ...
            'LineWidth', plotLineThickness, ...
            'Color', approach.colorAnimated);
        hLines{plotID}.HandleVisibility = 'off';

        % Add the markers
        hMarkers{plotID} =  plot(theta1Range, dataMarker, ...
            'LineStyle', 'none', ...
            'Marker', approach.marker, ...
            'MarkerSize', approach.markerSize, ...
            'Color', [0,0,0],... %approach.colorBW, ...
            'DisplayName', approachName);
        hMarkers{plotID}.HandleVisibility = 'off';

        % Dummy line with both line and markers for the legend
        hLDummy = plot(2, 2, ...
            'LineStyle', approach.lineStyle, ...
            'LineWidth', plotLineThickness, ...
            'Marker', approach.marker, ...
            'Color', approach.colorAnimated, ...
            'DisplayName', approachName);
        
 
    end

    % Set limits for axes
    xlim([min(theta1Range), max(theta1Range)])
    yLims = plotDescrStruct.yLims(currentParTable, paramName);
    ylim(yLims)

    % Add vertical line if plot is relative
    if plotDescrStruct.relative
        hV = yline(1);
        hV.HandleVisibility='off';
    end

    % Add a title and left-justify it
    switch plotID 
        case 1
            ttl = title('MSE relative to no averaging');
        case 2
            ttl = title('Bias');
        case 3
            ttl = title('Variance relative to no averaging');
    end 
    ttl.Units = 'Normalize';
    ttl.Position(1) = 0;
    ttl.HorizontalAlignment = 'left';

    % Add x-axis labels
    xlabel('Parameter space')
    ha(plotID).XTickLabel = ha(plotID).XTick;
     
    % Add y-ticks
    yTicks = determineTicks(yLims, 4, 10);
    ha(plotID).YTick = yTicks;
    ha(plotID).YTickLabel = yTicks;

    % Add legend based on the bias plot
    if plotID == 2
        legend(["No averaging", "Equal weights", "Optimal weights"]);
    end 

    % Add suptitle
    sgtitle('Optimal unit averaging leads to lower MSE and better bias-variance trade-off')
end

% Dynamically draw the fixed-N estimator
% Set the number of frames (drawing-pause-delete)
pauseFrames = gridMultiplier * 4;  
totalFrames = length(thetaGridMSE) + pauseFrames + length(thetaGridMSE);

% Initialize an array to store frames
frames = [];

% Animation loop
for i = 1:totalFrames
    if i <= length(thetaGridMSE)
        % Drawing phase

        % Extend lines
        set(hLines{1}, 'XData', thetaGridMSE(1:i), ...
            'YData', dataInterp{1}(1:i))
        set(hLines{2}, 'XData', thetaGridMSE(1:i), ...
            'YData', dataInterp{2}(1:i))
        set(hLines{3}, 'XData', thetaGridMSE(1:i), ...
            'YData', dataInterp{3}(1:i)) 

        % Handle adding markers when the line gets to them
        if rem(i-gridMultiplier, gridMultiplier) == 0
            markID = 1+(i-gridMultiplier)/gridMultiplier;
            set(hMarkers{1}, 'XData', theta1Range(1:markID), ...
                'YData',  currentLineData{1}(1:markID));
            set(hMarkers{2}, 'XData', theta1Range(1:markID), ...
                'YData',  currentLineData{2}(1:markID));
            set(hMarkers{3}, 'XData', theta1Range(1:markID), ...
                'YData',  currentLineData{3}(1:markID));
              
        end

    elseif i <= length(thetaGridMSE) + pauseFrames  
        % Do nothing, just hold the current state    

    else
        % Deletion phase

        % Rescale the current position being erased
        eraseIndex =  i-length(thetaGridMSE) - pauseFrames;
        
        % Clear lines
        set(hLines{1}, 'XData', thetaGridMSE(eraseIndex:end), ...
            'YData', dataInterp{1}(eraseIndex:end))
        set(hLines{2}, 'XData', thetaGridMSE(eraseIndex:end), ...
            'YData', dataInterp{2}(eraseIndex:end))
        set(hLines{3}, 'XData', thetaGridMSE(eraseIndex:end), ...
            'YData', dataInterp{3}(eraseIndex:end))
        
        % Clean up markers
        if rem(eraseIndex, gridMultiplier) == 0
            markID = (eraseIndex+gridMultiplier)/gridMultiplier;
            set(hMarkers{1}, 'XData', theta1Range(markID:end), ...
                'YData',  currentLineData{1}(markID:end));
            set(hMarkers{2}, 'XData', theta1Range(markID:end), ...
                'YData',  currentLineData{2}(markID:end));
            set(hMarkers{3}, 'XData', theta1Range(markID:end), ...
                'YData',  currentLineData{3}(markID:end));
              
        end

        % Clean up the plot on the last frame
        if i == totalFrames
            set(hLines{1}, 'XData', NaN, 'YData', NaN)
            set(hLines{2}, 'XData', NaN, 'YData', NaN)
            set(hLines{3}, 'XData', NaN, 'YData', NaN)
            set(hMarkers{1}, 'XData', NaN, 'YData', NaN)
            set(hMarkers{2}, 'XData', NaN, 'YData', NaN)
            set(hMarkers{3}, 'XData', NaN, 'YData', NaN)
        end
    end
    drawnow; 
    
    % Capture the frame
    frame = getframe(gcf);
    frames = [frames; frame];
end

drawSpeed = 0.034; % Speed of drawing the line 
eraseSpeed = 0.02;
% Write frames to GIF
filename = 'results/figures/animated_simplified_results.gif';
for k = 1:length(frames)
    im = frame2im(frames(k));
    [imind, cm] = rgb2ind(im, 256);
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, ...
            'DelayTime', drawSpeed);
    elseif k <= length(thetaGridMSE) + pauseFrames
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', ...
            'DelayTime', drawSpeed);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', ...
            'DelayTime', eraseSpeed);
    end
end


%% Finish export

close all

%% Auxiliary functions

function approachPos = findApproach(methodsArray, approachToFindShortName)
% findApproach Finds the index of an averaging approach by its short name.
%
% This function searches through a cell array of structs (`methodsArray`)
% to locate an averaging approach whose `shortName` field matches the
% specified string (`approachToFindShortName`). If found, it returns the
% index of the first match. If no match is found, it returns NaN.
%
% Args:
%     methodsArray (cell array of structs): 
%         Each struct represents an averaging approach and must have a 
%         `.shortName` field containing the short name of the approach.
%     approachToFindShortName (string): 
%         The short name of the approach to find.
%
% Returns:
%     approachPos (integer or NaN): 
%         The index of the first matching approach in `methodsArray`. 
%         Returns NaN if no match is found.

% Initialize the output to NaN
approachPos = NaN;

% Loop through each approach in the methods array
for approachIdx = 1:length(methodsArray)
    % Compare the shortName field with the target name
    if strcmp(methodsArray{approachIdx}.shortName, approachToFindShortName)
        % If a match is found, assign the index and exit the loop
        approachPos = approachIdx;
        break;
    end
end
end


function n = order(val, base)
% order Computes the order of magnitude of a number for a specified base.
%
% This function determines the order of magnitude of a given number `val` 
% with respect to a specified arithmetic base. If the base is not provided, 
% it defaults to 10. For example, the order of magnitude for 0.002 in base 
% 10 is -3, and for 1.3e6, it is 6.
%
% Args:
%     val (float): 
%         The number for which the order of magnitude is to be determined. 
%         Must be a non-zero value.
%     base (integer, optional): 
%         The base for the logarithmic computation. Defaults to 10 if not 
%         specified.
%
% Returns:
%     n (integer): 
%         The order of magnitude of `val` in the given `base`.

% Check if the base is provided, else default to base 10
if nargin < 2
    base = 10;
end

% Compute the order of magnitude using the logarithm and floor functions
n = floor(log(abs(val)) / log(base));
end


function ticks = determineTicks(yLims, minTicks, maxTicks)
% determineTicks Generates a vector of tick marks for axis scaling.
%
% This function calculates a suitable set of tick marks for a plot axis 
% based on given limits (`yLims`), while ensuring the number of ticks 
% falls within a specified range (`minTicks` to `maxTicks`). The ticks 
% are spaced logarithmically or linearly based on the scale of the input.
%
% Args:
%     yLims (vector of 2 floats): 
%         The axis limits [min, max], where `yLims(1)` is the lower limit 
%         and `yLims(2)` is the upper limit of the axis.
%     minTicks (integer, optional): 
%         Minimum number of ticks to generate. Defaults to 2.
%     maxTicks (integer, optional): 
%         Maximum number of ticks to generate. Defaults to 7.
%
% Returns:
%     ticks (vector): 
%         A vector of tick values, ordered from the lower to the upper axis . 
%         limit. The number of ticks is guaranteed to be between 
%         `minTicks` and `maxTicks`.

% Set default values for optional arguments
if nargin < 2
    minTicks = 2;
end
if nargin < 3
    maxTicks = 7;
end

% Determine the scale of the axis to calculate tick step size
scale = order(yLims(2)) - 1; 

% Round the lower limit to the nearest step based on scale
low = round(yLims(1), -scale);
high = yLims(2); % Upper limit remains unchanged

% Initialize step multiplier and generate initial tick vector
c = 1; 
ticks = low:(c * 10^scale):high;

% Adjust the step size iteratively to ensure the number of ticks is valid
while length(ticks) > maxTicks || length(ticks) < minTicks
    if length(ticks) > maxTicks
        % Reduce the number of ticks by increasing step size
        c = 2 * c;
    else
        % Increase the number of ticks by reducing step size
        c = 0.5 * c;
    end
    % Recompute the tick vector with the adjusted step size
    ticks = low:(c * 10^scale):high;
end
end


function plotName = ...
    approachPlotName(approach, approachRenames, approachShortName)
% approachPlotName Determines the display name for a given approach.
%
% This function returns the appropriate name to be used in plots for an 
% averaging approach. If the `approachShortName` is found in the 
% `approachRenames` struct, the corresponding renamed value is used. 
% Otherwise, the `longName` field from the `approach` struct is returned.
%
% Args:
%     approach (struct): 
%         A struct representing an approach, which must include a 
%         `longName` field  containing its default name.
%     approachRenames (struct): 
%         A struct where field names are short names of approaches, and  
%         their values are the custom names to be used in plots.
%     approachShortName (string): 
%         The short name of the approach for which the plot name is 
%         required.
% Returns:
%     plotName (string): 
%         The name to be displayed in the plot, which is either:
%         - The custom name from `approachRenames` (if available), or
%         - The `longName` field from the `approach` struct.
%

% Check if a custom name is defined in `approachRenames` for the short name
if isfield(approachRenames, approachShortName)
    % Use the custom name from approachRenames
    plotName = approachRenames.(approachShortName);
else
    % Fallback to the default long name from the approach struct
    plotName = approach.longName;
end
end


function plotArray = changeLinesToPlot(plotArray, linesToPlot)
% changeLinesToPlot Updates the plotting configuration with a new set of
% unit averaging approaches.
%
% This function modifies the `.approachesToPlot` field in each element of 
% the `plotArray` to include the provided `linesToPlot`. This is useful for
% dynamically changing the lines or methods that appear in plots.
%
% Args:
%     plotArray (cell array): A cell array of structs, where each struct 
%         represents a plot description and includes an `.approachesToPlot`
%         field.
%     linesToPlot (cell array of strings): A vector of strings specifying 
%         the short names of the unit averaging approaches to include in 
%         the plots.
%
% Returns:
%     plotArray (cell array): The updated plot descriptions with modified 
%         `.approachesToPlot` fields reflecting the new `linesToPlot`
%         configuration.

% Iterate through each plot configuration in the input array
for plotID = 1:length(plotArray)
    % Update the `approachesToPlot` field with the provided `linesToPlot`
    plotArray{plotID}.approachesToPlot = linesToPlot;
end
end


function figureSavingName = makeFigureFileName(...
    leadingString, simulationSetting, ...
    plotSavingName, paramName, ...
    valuesN, valuesT)
    % makeFigureFileName Generates a standardized figure name for exporting
    % figure
    %
    % Args:
    %     leadingString (string): String to insert at the beginning of the
    %         file name.
    %     simulationSetting (string): Name of the data-generating process.
    %     plotSavingName (string): Name characterizing figure contents.
    %     paramName (string): Name of parameter being plotted.
    %     valuesN (vector): Vector of cross-sectional sizes used.
    %     valuesT (vector): Vector of time series sizes used.
    %
    % Returns:
    %     figureSavingName (str): Constructed file name including all  
    %         parameters, saved in the 'results/simulation' folder. Format:
    %         'results/figures/[leadingString][plotSavingName]_...
    %         [simulationSetting]_[paramName]_N-[valuesN]_T-[valuesT]'
    
    % Numbers of dimensions of samples drawn
    numT = length(valuesT);
    numN = length(valuesN);

    % Generate parts of the title that depend on (N, T)
    titleN = "";
    for nID=1:numN
        titleN = titleN + "-" + valuesN(nID) ;
    end

    titleT = "";
    for tID=1:numT
        titleT = titleT + "-" + valuesT(tID);
    end
    
    % Create figure title
    figureSavingName = sprintf('results/figures/%s%s_%s_%s_N%s_T%s', ...
        leadingString, ...          % Insert the leading string
        plotSavingName, ...         % Plot name
        simulationSetting, ...      % Simulation context
        paramName, ...              % Parameter name
        titleN, titleT);            % Sample sizes and cross-sections 
end