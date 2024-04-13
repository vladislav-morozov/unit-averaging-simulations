%% Plotting parameters

% Specify approaches to appear in the main text
approachesToPlotMSEPaper = ["ind", "mg", "aic", ...
    "unrestr", "oracle_similar_10", "stein", "top_10_pct","cluster_coef_4"];

% Specify any renames required from the default values
approachPaperRenames.top10 = "10 most similar";



% Specify which approaches to plot for each section
% approachesToPlotMSE = ["ind", "mg", "aic", "mma", "unrestr",...
%     "cluster_coef_2", "cluster_coef_4", "top", "stein", "random10pct", "oracleSimilarity" ];

approachesToPlotFirstWeight = ["unrestr", "top", "oracleSimilarity", "stein"];

% Plotting parameters
plotLineThickness= 1.3;
plotWLines = 1000; % width in pixels
plotHLines = 1600; % height in pixels

% Grid of points for line plots. Markers will be placed at positions 
% theta1Range only
thetaGridMSE = min(theta1Range):0.01:max(theta1Range);

% Axis plotting limits
yLimsRelMSE.lambda = [0.5, 1.3];
yLimsRelMSE.beta = [0.9, 1.2];
yLimsRelMSE.longRun = [0.7, 1.2];
yLimsRelMSE.forecast = [0.8, 1.4];

% Set text interpreters to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex'); 


%% Line plots template 

approachesToPlotDummy = ["unrestr"];

% Relative MSE
linePlots{1}.data = mseTablesN;
linePlots{1}.approachesToPlot = approachesToPlotDummy;
linePlots{1}.firstColumnPlot = 2;
linePlots{1}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx}./dataTable{:, 1};
linePlots{1}.yLims = @(dataTable, paramName) ...
    [yLimsRelMSE.(paramName)];

linePlots{1}.relative = true;
linePlots{1}.suptitle = @(plotDescr) {'Averaging estimators, ' + ...
    plotDescr,  'Ratio of MSE to individual estimator'};
linePlots{1}.plotSavingName = 'relative_mse';

% Absolute MSE
linePlots{2}.data = mseTablesN;
linePlots{2}.approachesToPlot = approachesToPlotDummy;
linePlots{2}.firstColumnPlot = 1;
linePlots{2}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{2}.yLims = @(dataTable, paramName) [0, 1.5*max(dataTable.unrestr)];
linePlots{2}.relative = false;
linePlots{2}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', MSE';
linePlots{2}.plotSavingName = 'absolute_mse';

% Bias
linePlots{3}.data = biasTablesN;
linePlots{3}.approachesToPlot = approachesToPlotDummy;
linePlots{3}.firstColumnPlot = 1;
linePlots{3}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx});
%     @(dataTable, colIdx) (dataTable{:, colIdx}.^2)./(dataTable{:, 1}.^2);
linePlots{3}.yLims = @(dataTable, paramName) [1.5*min(dataTable.unrestr), ...
    1.5*max(dataTable.unrestr)];
%     1.5*max(  ((dataTable.unrestr).^2)./((dataTable{:, 1}.^2)))];
linePlots{3}.relative = false;
linePlots{3}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', bias';
linePlots{3}.plotSavingName = 'bias';

% Relative variance
linePlots{4}.data = varTablesN;
linePlots{4}.approachesToPlot = approachesToPlotDummy;
linePlots{4}.firstColumnPlot = 2;
linePlots{4}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx})./(dataTable{:, 1});
linePlots{4}.yLims = @(dataTable, paramName) [0, ...
    1.5*max(  (dataTable.unrestr./(dataTable{:, 1})))];
linePlots{4}.relative = true;
linePlots{4}.suptitle = @(plotDescr) 'Averaging estimators, ' + ...
    plotDescr + ...
    ', relative variance';
linePlots{4}.plotSavingName = 'relative_variance';


% Average weight of first unit
linePlots{5}.data = firstWeightNT;
linePlots{5}.approachesToPlot = approachesToPlotDummy;
linePlots{5}.firstColumnPlot = 1;
linePlots{5}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{5}.yLims = @(dataTable, paramName) [0, 1];
linePlots{5}.relative = false;
linePlots{5}.suptitle = @(plotDescr) {'Minimum MSE estimators, ' + ...
    plotDescr, ...
    'Average weight of target unit'};
linePlots{5}.plotSavingName = 'average_first_weight';


%% Main Text 1 and 2
% This sections generates the plots of the first two figures of section 4
% Fig 1. Relative MSE (lambda; all N; T=30,60, 600)
% Fig 2. Average own weight (lambda; N=50,150; T = 60,600)

% Set line plots to use lines for the paper
linePlots = changeLinesToPlot(linePlots, approachesToPlotMSEPaper);

% Parameter to report in the main text
parID = 1;

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
            axes(ha(plotNum))  
            
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


%% Main Text 3: Bias and Variance plots
% Plot bias and variance for only one combination of (N, T)
% Chosen values: N=150, T=60

% Set value of (N, T) to use
nID = 2; % N=150
tID = 2; % T=60
% Parameter: lambda
parID = 1;
paramName = paramArray{parID}.saveName;

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
    [.04 .04],[.07 0.08],[.05 .01]);

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
        approachName = approach.longName;
        
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


%% Full-color plots
% Plots for section Extra MSE results
% to be replaced by the final version of 

% Loop through plot descriptions
% for plotID = 1:1
for plotID = 1:length(linePlots)
    % Extract description of the plot
    plotDescrStruct = linePlots{plotID};
    
    % Extract the approaches for this plot
    approachesToPlot = plotDescrStruct.approachesToPlot;
    
    % Loop over parameter   s
    for parID = 1:numParams
        %     for parID = 1:1
        
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
                    approachName = approach.longName;
                    
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
                        'MarkerSize', 6.7, ...
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
                yLims = plotDescrStruct.yLims(currentParTable, paramName);
                ylim(yLims)
                
                % Add vertical line if appropriate
                if plotDescrStruct.relative
                    hV = yline(1);
                    hV.HandleVisibility='off';
                end
                
                % Add a title and left-justify it
                ttl = title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)));
                ttl.Units = 'Normalize';
                ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
                ttl.HorizontalAlignment = 'left';
                
                % x axis labels only on the bottom
                if tID == numT
                    xlabel('\lambda_1')
                    ha(plotNum).XTickLabel = ha(plotNum).XTick;
                end
                % y axis labels only on the left
                % y axis labels only on the left
                if nID == 1
                    yTicks = determineTicks(yLims);
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
        legend('Position', [0.85, 0.91, 0.1, 0.07]);
        
        % Save figure
        figureSavingName =  figureFolderName + "/" + ...
            plotDescrStruct.plotSavingName + "_" + ...
            simulationSetting+'_'+paramName +...
            "_N" + titleN + ...
            "_T" + titleT ;
        
        
        % Save both EPS and JPG versions
        %         saveas(p,figureSavingName,'epsc');
        %          saveas(p,figureSavingName,'jpg')
        set(gcf, 'PaperPosition', [0 0 10 12])    % can be bigger than screen
        print(gcf, figureSavingName, '-dpng', '-r300' );
    end
end






%% 
determineTicks(yLims, 3, 10)
%% 1D Plot of Weights

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
    
    % Loop over N and T
    for tID = 1:numT
        for nID = 1:numN
            % Extract weights
            weightsNW = weightsTablesNT{nID, tID}.(paramName).(approachWeightPlot);
            
            % Create subplot
            ax = subplot(numT, numN, (tID-1)*numN+nID);
            hold(ax, 'on');
            % Plot as individual lines
            for targetID = 1:2:length(theta1Range)
                lineName = "\lambda_1 = " + num2str(theta1Range(targetID));
                plot(thetaGrid, weightsNW(targetID, :),...
                    'DisplayName', lineName, ...
                    'color', [targetID/length(theta1Range), 0, 0])
            end
            
            if tID == numT && nID == numN
                legend()
            end
            
            
            title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)))
            xlabel('\lambda_{alt}')
        end
    end
    suptitle("Average weight of unit with \lambda_{alt} for estimating parameter of unit 1, " + approachLongName + ", " + paramArray{parID}.plotDescr)
    % Save figure
    figureSavingName =  figureFolderName + "/" + ...
        "weights_1d_" + ...
        simulationSetting+'_'+paramName + ...
        "_N" + titleN + ...
        "_T" + titleT ;
    
    % Save both EPS and JPG versions
    saveas(p,figureSavingName,'epsc');
    saveas(p,figureSavingName,'jpg')
end

%% 2D Plot of Weights
% Use the 'hot' colormap


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
    
    % Loop over N and T
    for tID = 1:numT
        for nID = 1:numN
            % Extract weights
            weightsNW = weightsTablesNT{nID, tID}.(paramName).(approachWeightPlot);
            
            % Create subplot
            ax = subplot(numT, numN, (tID-1)*numN+nID);
            hold(ax, 'on');
            surf(X, Y, weightsNW)
            
            
            
            
            
            title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)))
            ylabel('\lambda_1')
            xlabel('\lambda_{alt}')
        end
    end
    suptitle("Average weight of unit with \lambda_{alt} for estimating parameter of unit 1, " + approachLongName + ", " + paramArray{parID}.plotDescr)
    % Save figure
    figureSavingName =  figureFolderName + "/" + ...
        "weights_2d_" + ...
        simulationSetting+'_'+paramName + ...
        "_N" + titleN + ...
        "_T" + titleT ;
    
    % Save both EPS and JPG versions
    saveas(p,figureSavingName,'epsc');
    saveas(p,figureSavingName,'jpg')
end


 
%% Test space


changeLinesToPlot(linePlots, approachesToPlotMSEPaper) 


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
% Author Ivar Smith
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