% ===========================================================
% File: exportFigures.m
% Description: This script exports the figures based on simulations for the
%              mean squared error, bias, variance, weights, and
%              unrestricted units of unit averaging
%
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
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

% Load line plots
chooseLinePlots

%% Main Text Plots
% This sections generates the plots of the first two figures of section 4
% Fig 1. Relative MSE (lambda; all N; T=30,60, 600)
% Fig 2. Example bias and variance for  single value of (N, T)
% Fig 3. Average own weight (lambda; N=50,150; T = 60,600)

% Set line plots to use lines for the paper
linePlots = changeLinesToPlot(linePlots, approachesToPlotMSEPaper);

% Parameter to report in the main text
parID = 1;

% Block 1: figures 1 and 3
% Plot matrices to generate for the main text
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
textPlotUpperMargin{1} = 0.02;
textPlotUpperMargin{2} = 0.16;

% Create plots
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
            'Position', [50 50 1900 1200]);
    else
        figure('visible','off',...
            'Renderer', 'painters',...
            'Position', [50 50 plotWLines plotHLines]);
    end
    
    % Create tight axes using above parameters
    [ha, ~] = tight_subplot(numTPlot,numNPlot, ...
        [.052 .023],[.07 textPlotUpperMargin{textPlotID}],[.05 .01]);
    
    % Loop over N and T included in this plot
    tPlotID = 0;
    nPlotID = 0;
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
            
            % x axis labels only on the bottom
            if tID == numT
                xlabel('$\lambda_1$')
                ha(plotNum).XTickLabel = ha(plotNum).XTick;
            end
            % y axis labels only on the left
            if nID == 1
                yTicks = determineTicks(yLims, 4, 10);
                ha(plotNum).YTick = yTicks;
                ha(plotNum).YTickLabel = yTicks;
            end
            % Enforce same ticks as in first column
            ha(plotNum).YTick = yTicks;
        end
    end
    
    % Include suptitle
    suptitle(plotDescrStruct.suptitle(paramArray{parID}.plotDescr))
    
    % Add legend
    legend('Position', [0.83, 0.88, 0.1, 0.06]);
    
    % Name for the saved figure
    figureSavingName =  figureFolderName + "/BW" + ...
        plotDescrStruct.plotSavingName + "_" + ...
        simulationSetting+'_'+paramName +...
        "_N" + titleN + ...
        "_T" + titleT ;
    
    % Save both EPS and PNG versions
    set(gcf, 'PaperPosition', textPlotSavingSize{textPlotID})
    print(gcf, figureSavingName, '-depsc' );
    print(gcf, figureSavingName, '-dpng', '-r300' );
end

%% Experimental BW bias and variance plot
% Parameter to report in the main text

nIDs = [2 2, 3, 3];
tIDs = [1 2 1 2 ];

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
[ha, pos] = tight_subplot(2, 2, ...
    [.04 .04],[.07 0.08],[.05 .01]);

% plotID 3 and 4 correspond to bias and variance
for plotNum = 1:4
    nID = nIDs(plotNum); %
    tID = tIDs(plotNum); % T=60
    % Axes ID to use
    plotID = 3+(mod(plotNum, 2)==0);
    axes(ha(plotNum)) %, 'Parent', p);
    
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

% Name for the saved figure
figureSavingName =  figureFolderName + "/BW" + ...
    "bias_variance_2" + "_" + ...
    simulationSetting+'_'+paramName +...
    "_N_150" + ...
    "_T_60" ;

% Save both EPS and JPG versions
set(gcf, 'PaperPosition', [0, 0, 13, 4]*0.8)
print(gcf, figureSavingName, '-depsc' );
print(gcf, figureSavingName, '-dpng', '-r300' );


%% Block 2: figure 2
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
[ha, pos] = tight_subplot(1, 2, ...
    [.04 .04],[.14 0.08],[.05 .01]);

% plotID 3 and 4 correspond to bias and variance
for plotID = 3:4
    % Axes ID to use
    plotNum = plotID-2;
    axes(ha(plotNum)) %, 'Parent', p);
    
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

% Name for the saved figure
figureSavingName =  figureFolderName + "/BW" + ...
    "bias_variance" + "_" + ...
    simulationSetting+'_'+paramName +...
    "_N_150" + ...
    "_T_60" ;

% Save both EPS and JPG versions
set(gcf, 'PaperPosition', [0, 0, 13, 4]*0.8)
print(gcf, figureSavingName, '-depsc' );
print(gcf, figureSavingName, '-dpng', '-r300' );


%% Full-color line plots for the Online Appendix
% This sections generates the plots of the Online Appendix
% 1. MSE, own average weights for all the approaches considered in the
%    paper, extra parameters, including T=180. In color.
% 2. MSE for groups of tuning parameter choices



% for lineSetID = 1:length(lineSets)
for lineSetID = 9:9
    
    % Set line plots to use
    linePlots = changeLinesToPlot(linePlots, lineSets{lineSetID});
    
    % Loop through plot descriptions
    for plotID = 1:1
%     for plotID = 1:length(linePlots)
        % Extract description of the plot
        plotDescrStruct = linePlots{plotID};
        
        % Extract the approaches for this plot
        approachesToPlot = plotDescrStruct.approachesToPlot;
        
        % Loop over parameter   s
%         for parID = 1:numParams
        for parID = 1:1
            
            % For main set of lines, generate MSE and weights for all, and
            % bias and variance only for lambda
            if parID > 1 && ismember(plotID, [2, 3, 4])% && lineSetID == 1
                continue
            end
            % For comparison lines, only generate the MSE
            %             if lineSetID > 1 && plotID>1
            %                 continue
            %             end
            
            
            
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
                    
                    % Patch bug: in the  implementation of the large-N top
                    % method
                    currentParTable.top_50_pct = ...
                        currentParTable.unrestr;
                    if nID==1 && tID == 4
                        currentParTable.top_25 = ...
                            currentParTable.unrestr;
                    end
                    
                    % Extract approaches present in this table
                    approachesPresent = currentParTable.Properties.VariableNames;
                    
                    % Extract number of approaches
                    numApproaches = size(currentParTable, 2);
                    
                    % Plot
                    hold(ha(plotNum), 'on');
                    for approachID = plotDescrStruct.firstColumnPlot:numApproaches
                        
                        % Get approach short name
                        approachShortName = string(approachesPresent{approachID});
                        
                        % Check if approach is present in the list of
                        % approaches to plot
                        approachPresent = ismember(approachShortName ,approachesToPlot);
                        % Skip if this approach was not available
                        if ~approachPresent
                            continue
                        end
                        
                        % Extract information
                        approachCoord = findApproach(allMethodsArray, approachShortName);
                        approach = allMethodsArray{approachCoord};
                        approachName = approachPlotName(approach, ...
                            nameSets{lineSetID}, ...
                            approachShortName);
                        
                        
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
                        
                        % Dummy line
                        hLDummy = plot(2, 2, ...
                            'LineStyle', approach.lineStyle, ...
                            'LineWidth', plotLineThickness, ...
                            'Marker', approach.marker, ...
                            'Color', approach.color, ...
                            'DisplayName', approachName);
                    end
                    %                 hold(ax,'off');
                    % Set limits for axes
                    xlim([min(theta1Range), max(theta1Range)])
                    if nID == 1
                        yLims = plotDescrStruct.yLims(currentParTable, paramName);
                    end
                    ylim(yLims)
                    
                    % Add vertical line if appropriate
                    if plotDescrStruct.relative
                        hV = yline(1);
                        hV.HandleVisibility='off';
                    end
                    if plotID == 7
                       hV = yline(0);
                        hV.HandleVisibility='off'; 
                    end
                    
                    % Add a title and left-justify it
                    ttl = title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)));
                    ttl.Units = 'Normalize';
                    ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
                    ttl.HorizontalAlignment = 'left';
                    
                    % x axis labels only on the bottom
                    if tID == numT
                        xlabel('$\lambda_1$')
                        ha(plotNum).XTickLabel = ha(plotNum).XTick;
                    end
                    % y axis labels only on the left
                    if nID == 1
%                         if min(yLims)>=0
                        yTicks = determineTicks(yLims, 3, 10);
%                         else
%                          yTicks = -determineTicks(-yLims, 3, 10);   
%                         end
                        ha(plotNum).YTick = yTicks;
                        ha(plotNum).YTickLabel = yTicks;
                    end
                    % Enforce same ticks as in first column
                    ha(plotNum).YTick = yTicks;
                end
            end
            
            
            % Include suptitle
            sgtitle(plotDescrStruct.suptitle(paramArray{parID}.plotDescr))
            
            % Add legend
            legend('Position', [0.845, 0.9, 0.08, 0.066]);
            
            % Save figure
            figureSavingName =  figureFolderName + "/" + ...
                plotName{lineSetID} + "_" +...
                plotDescrStruct.plotSavingName + "_" + ...
                simulationSetting+'_'+paramName;
            %                 "_N" + titleN + ...
            %                 "_T" + titleT ;
            
            
            % Save both EPS and JPG versions
            %         saveas(p,figureSavingName,'epsc');
            %          saveas(p,figureSavingName,'jpg')
            % Save both EPS and PNG versions
            set(gcf, 'PaperPosition', [0, 0, 10, 12])
            print(gcf, figureSavingName, '-depsc' );
            print(gcf, figureSavingName, '-dpng', '-r300' );
        end
    end
end


%% 1D Plot of Weights 

yLims = [0, 0.1];

% Approach to plot
approachWeightPlot = "unrestr";
for approachID = 1:length(allMethodsArray)
    if allMethodsArray{approachID}.shortName == approachWeightPlot
        break
    end
end

approach = allMethodsArray{approachID};
approachLongName = approach.longName;
for parID = 1:numParams
% for parID = 1:1
    % Parameter name
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
                % Create plotting grid
            [ha, ~] = tight_subplot(numT,numN,[.04 .023],[.07 .01],[.05 .01]);
            
    % Loop over N and T
    for tID = 1:numT
        for nID = 1:numN
            % Extract weights
            weightsNW = weightsTablesNT{nID, tID}.(paramName).(approachWeightPlot);
            
            % Create subplot
            % Create subplot
            plotNum = (tID-1)*numN+nID;
            axes(ha(plotNum))
            
            % Plot as individual lines
            for targetID = 1:3:length(theta1Range)
                lineName = "$\lambda_1$ = " + num2str(theta1Range(targetID));
                plot(thetaGrid, weightsNW(targetID, :),...
                    'DisplayName', lineName, ...
                    'LineWidth', 1+targetID/(2*length(theta1Range)), ...
                    'color', [targetID/length(theta1Range), 0, 0])
                hold(ha(plotNum), 'on');
                xlim([min(theta1Range), max(theta1Range)]);
                %                     'LineWidth', 2-targetID/(2*length(theta1Range)), ...
            end
            
            if tID == numT && nID == numN
                legend()
            end
            ylim(yLims)
            % x axis labels only on the bottom
            if tID == numT
                xlabel('$\lambda_1$')
                ha(plotNum).XTickLabel = ha(plotNum).XTick;
            else
                ha(plotNum).XTickLabel = [];
            end
            % y axis labels only on the left
            if nID == 1
                yTicks = determineTicks(yLims, 4, 10);
                ha(plotNum).YTick = yTicks;
                ha(plotNum).YTickLabel = yTicks;
            else
                ha(plotNum).YTickLabel = [];
            end
            % Enforce same ticks as in first column
            ha(plotNum).YTick = yTicks;
            
            % Add a title and left-justify it
            ttl = title("N="+num2str(valuesN(nID))+...
                ", T="+num2str(valuesT(tID)));
            ttl.Units = 'Normalize';
            ttl.Position(1) = 0;
            ttl.HorizontalAlignment = 'left';
            
            if nID == 1
                ylabel('$\lambda_1$')
            end
            if tID == numT
                xlabel('$\lambda_{alt}$')
            end
        end
    end
    suptitle({"Average weight of unit with $\lambda_{alt}$ for estimating parameter of unit 1",...
        approachLongName + ", " + paramArray{parID}.plotDescr})
    % Save figure
    figureSavingName =  figureFolderName + "/" + ...
        "weights_1d_" + ...
        simulationSetting+'_'+paramName;
    
            set(gcf, 'PaperPosition', [0, 0, 10, 12])
            print(gcf, figureSavingName, '-depsc' );
            print(gcf, figureSavingName, '-dpng', '-r300' );
end

%% 2D Plot of Weights and Unrestricted Units


gridPlots{1}.data = weightsTablesNT;
% gridPlots{1}.approachesToPlot = approachesToPlotMSEPaper([4, 5, 7, 8]);
gridPlots{1}.approachesToPlot = approachesAlt;
gridPlots{1}.logScale = true;
gridPlots{1}.cmapN1 = 150;
gridPlots{1}.cmapN2 = 45;
gridPlots{1}.cmapFlip = true;
gridPlots{1}.caxisLims = [0.001, 1];
gridPlots{1}.plotSavingName = "weights_2d";
gridPlots{1}.suptitle = ...
    "Average weight of a unit with $\lambda_{alt}$ in estimating $\lambda_1$";

gridPlots{2}.data = unitsUnrestrNT;
% gridPlots{2}.approachesToPlot = approachesToPlotMSEPaper([4, 5, 7, 8]);
gridPlots{2}.approachesToPlot = approachesAlt;
% gridPlots{2}.logScale = true;
% gridPlots{2}.cmapN1 = 25;
% gridPlots{2}.cmapN2 = 120;
gridPlots{2}.logScale = false;
gridPlots{2}.cmapN1 = 150;
gridPlots{2}.cmapN2 = 20;
gridPlots{2}.cmapFlip = true;
gridPlots{2}.caxisLims = [0, 1];
gridPlots{2}.plotSavingName = "unrestr_2d";
gridPlots{2}.suptitle = ...
    "Probability of a unit with $\lambda_{alt}$ being unrestricted in estimating $\lambda_1$";


approachesToPlotGridWeights = approachesToPlotMSEPaper([4, 5, 7, 8]);

for plotID = 1:length(gridPlots)
    % for plotID = 1:1
    for approachID = 1:length(gridPlots{plotID}.approachesToPlot)
        %     for approachID = 2:4
        approachWeightPlot = gridPlots{plotID}.approachesToPlot(approachID);
        for candApproachID = 1:length(allMethodsArray)
            if allMethodsArray{candApproachID}.shortName == approachWeightPlot
                break
            end
        end
        
        
        approach = allMethodsArray{candApproachID};
        approachName = approachPlotName(approach, ...
            approachGridRenames, ...
            approachWeightPlot);
%         for parID = 1:numParams
                    for parID = 1:1
            % Parameter name
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
            [ha, ~] = tight_subplot(numT,numN,[.04 .023],[.07 .01],[.05 .01]);
            
            % Loop over N and T
            for tID = 1:numT
                for nID = 1:numN
                    % Extract weights
                    dataNW = gridPlots{plotID}.data{nID, tID}.(paramName).(approachWeightPlot);
                    
                    % Create subplot
                    plotNum = (tID-1)*numN+nID;
                    axes(ha(plotNum))
                    
                    pcolor(X, Y, dataNW)
                    %             set(ax,'zscale','log')
                    
                    % White for smallest values, blue for intermediate, orange for
                    % largest
                    colormap(...
                        orangeBlueMaps(...
                        gridPlots{plotID}.cmapN1,...
                        gridPlots{plotID}.cmapN2,...
                        gridPlots{plotID}.cmapFlip))
                    if nID == 1
                        zLimits = get(gca,'ZLim') ;
                    end
                    caxis(gridPlots{plotID}.caxisLims)
                    tmpAspect=daspect(); % get the aspect ratio of the axes scales
                    daspect(tmpAspect([1 2 2]));
                    
                    %             zlim(zLimits);
                    if gridPlots{plotID}.logScale
                        set(gca,'ColorScale','log')
                    end
                    if nID == numN
                        colorbar
                        
                    end
                    
                    
                    
                    
                    % Add a title and left-justify it
                    ttl = title("N="+num2str(valuesN(nID))+...
                        ", T="+num2str(valuesT(tID)));
                    ttl.Units = 'Normalize';
                    ttl.Position(1) = 0;
                    ttl.HorizontalAlignment = 'left';
                    
                    if nID == 1
                        ylabel('$\lambda_1$')
                    end
                    if tID == numT
                        xlabel('$\lambda_{alt}$')
                    end
                end
            end
            suptitle({gridPlots{plotID}.suptitle,
                approachName + ", " + paramArray{parID}.plotDescr})
            % Save figure
            figureSavingName =  figureFolderName + "/" + ...
                gridPlots{plotID}.plotSavingName+ '_' + ...
                simulationSetting+'_'+...
                approachWeightPlot + '_' +paramName;
            
            set(gcf, 'PaperPosition', [0, 0, 10, 12])
            print(gcf, figureSavingName, '-depsc' );
            print(gcf, figureSavingName, '-dpng', '-r300' );
            
        end
    end
end

%% Functions used for plotting in this script
% Functions used
% 1. findApproach -- finding coordinate of approach using its name
% 2. order -- computes order of number
% 3. determineTicks -- determines ticks for plots
% 4. approachPlotName -- returns name for plotting


function approachPos = findApproach(methodsArray, approachToFindShortName)
% findApproach Finds the cooordinate of an averaging approach using its
% short name
% Args:
%       1. methodsArray -- cell array of structs. Each struct is an
%       averaging approach with a .shortName field
%       2. approachToFindShortName -- string, short name of the approach to
%       find.
%
% Outputs:
%       1. approachPos -- integer or NaN. If approach in the methodsArray,
%       then returns the index of the (first matching) approach. If
%       approach not in the array, returns NaN.

% Default output: NaN
approachPos = NaN;

% Loop through approaches
for approachIdx = 1:length(methodsArray)
    % If shortName matches the target, set the output and finish the search
    if methodsArray{approachIdx}.shortName == approachToFindShortName
        approachPos = approachIdx;
        break
    end
    
end
end


function n = order(val, base)
% order Order of magnitude of number for specified base. Default base is 10
% order(0.002) will return -3., order(1.3e6) will return 6. 
if nargin < 2
    base = 10;
end
n = floor(log(abs(val))./log(base));
end


function ticks = determineTicks(yLims, minTicks, maxTicks)
% determineTicks Determines the ticks for the plots of this script.
%
% Args:
%       1. yLims -- 2 vector of limits for the axis
%       2. minTicks -- integer, optional, default value 2. Minimum number
%                      of ticks
%       3. maxTicks -- integer, optional, default value 7. Maximum number
%                      of ticks
% Outputs:
%       1. ticks -- a vector of length between minTicks and maxTicks.
%                   Contains ticks ordered tick values

% Set default values
if nargin<2
    minTicks = 2;
    maxTicks = 7;
end

% Extract scale of axis to creat steps
scale = order(yLims(2))-1;

% Set endpoints for ticks
low = round(yLims(1),-(scale));
high = yLims(2);

% Try setting the number of ticks by default
c = 1;
ticks = low:(c*10^scale):high;

% Increase or decrease the number of ticks until producing a valit output
while length(ticks)> maxTicks || length(ticks)<minTicks
    % Too long: increase coarseness
    if length(ticks)> maxTicks
        c = 2*c;
        % Too short: use finer steps
    else
        c = 0.5*c;
    end
    % Compute ticks with new step
    ticks = low:(c*10^scale):high;
end
end


function plotName = approachPlotName(...
    approach, ...
    approachRenames, ...
    approachShortName)
% approachPlotName Returns the plotting name for a given approach. If
% approach is in approachRenames, the renamed name will be used. Otherwise
% defauts to the longName in approach
%
% Args:
%   1. approach -- struct with field .longName
%   2. approachRenames -- struct. Field names are approach short names,
%      values are names to display
%   3. approachShortName -- string
%
% Outputs:
%   1. plotName -- either name in approachRenames or the longName from
%                  approach

% If approach is to be renamed, use the new name
if isfield(approachRenames, approachShortName)
    plotName = approachRenames.(approachShortName);
else
    plotName = approach.longName;
end
end


function linePlots = changeLinesToPlot(linePlots, linesToPlot)
% changeLinesToPlot Replaces the .approachesToPlot field in every member of
% linePlots with linesToPlot
for plotID = 1:length(linePlots)
    linePlots{plotID}.approachesToPlot = linesToPlot;
end
end