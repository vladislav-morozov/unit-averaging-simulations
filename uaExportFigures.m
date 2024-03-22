
% Plotting parameters
plotLineThickness= 1.1;
plotW = 1600; % width in pixels
plotH = 400; % height in pixels



%% Line plots

% Relative MSE
linePlots{1}.data = mseTablesN;
linePlots{1}.firstColumnPlot = 2;
linePlots{1}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx}./dataTable{:, 1};
linePlots{1}.yLims = @(dataTable) [yMinRelMSEPlot, yMaxRelMSEPlot];
linePlots{1}.relative = true;
linePlots{1}.suptitle = @(plotDescr) 'Averaging estimator, ' + ...
                plotDescr + ...
               ', ratio of MSE to individual estimator';
linePlots{1}.plotSavingName = 'relative_mse';   

% Absolute MSE
linePlots{2}.data = mseTablesN;
linePlots{2}.firstColumnPlot = 1;
linePlots{2}.dataTransform = ...
    @(dataTable, colIdx) dataTable{:, colIdx};
linePlots{2}.yLims = @(dataTable) [0, 1.5*max(dataTable.unrestr)];
linePlots{2}.relative = false;
linePlots{2}.suptitle = @(plotDescr) 'Averaging estimator, ' + ...
                plotDescr + ...
                ', MSE';
linePlots{2}.plotSavingName = 'absolute_mse';

% Relative square bias
linePlots{3}.data = biasTablesN;
linePlots{3}.firstColumnPlot = 2;
linePlots{3}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx}.^2)./(dataTable{:, 1}.^2);
linePlots{3}.yLims = @(dataTable) [0, ...
    1.5*max(  ((dataTable.unrestr).^2)./((dataTable{:, 1}.^2)))];
linePlots{3}.relative = true;
linePlots{3}.suptitle = @(plotDescr) 'Averaging estimator, ' + ...
                plotDescr + ...
                ', relative squared bias';
linePlots{3}.plotSavingName = 'relative_sq_bias';

% Relative variance
linePlots{4}.data = varTablesN;
linePlots{4}.firstColumnPlot = 2;
linePlots{4}.dataTransform = ...
    @(dataTable, colIdx) (dataTable{:, colIdx})./(dataTable{:, 1});
linePlots{4}.yLims = @(dataTable) [0, ...
    1.5*max(  (dataTable.unrestr./(dataTable{:, 1})))];
linePlots{4}.relative = true;
linePlots{4}.suptitle = @(plotDescr) 'Averaging estimator, ' + ...
                plotDescr + ...
                ', relative variance';
linePlots{4}.plotSavingName = 'relative_variance';


%%
 

% Loop through plot descriptions
for plotID = 1:length(linePlots) 
    % Extract description of the plot
    plotDescrStruct = linePlots{plotID};
    
    % Loop over parameters
    for parID = 1:numParams
       
        % Parameter name
        paramName = paramArray{parID}.saveName;
       
        % Create figure
        if plotQuietly ~= 1
            p = figure('Renderer', 'painters', ...
                'Position', [50 50 plotW plotH]);
        else
            p = figure('visible','off',...
                'Renderer', 'painters',...
                'Position', [50 50 plotW plotH]); 
        end
        
        % Loop over N and T
        for tID = 1:numT
            for nID = 1:numN
                % Create subplot
                ax = subplot(numT, numN, (tID-1)*numN+nID);
                
                % Extract data
                currentParTable = plotDescrStruct.data{nID, tID}.(paramName);
                
                % Extract number of approaches
                numApproaches = size(currentParTable, 2); 
                
                % Plot
                hold(ax, 'on');
                for approachID = plotDescrStruct.firstColumnPlot:numApproaches
                    %                 if ~plotMethod(approachID)
                    %                    continue
                    %                 end
                    %
                    approach = methodsForPlotting{approachID};
                    approachName = approach.longName;
                    
                    currentLineData = ...
                        plotDescrStruct.dataTransform(...
                        currentParTable, approachID);

                    plot(theta1Range, currentLineData, ...
                        'LineStyle', approach.lineStyle, ...
                        'Marker', approach.marker, ...
                        'LineWidth', plotLineThickness, ...
                        'Color', approach.color, ...
                        'DisplayName', approachName)
                    
                    
                end
                hold(ax,'off');
                % Set limits for axes
                xlim([min(theta1Range), max(theta1Range)])
                ylim(plotDescrStruct.yLims(currentParTable))
                
                % Add vertical line if appropriate
                if plotDescrStruct.relative
                    h = yline(1);
                    h.HandleVisibility='off';
                end
                
                title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)))
                xlabel('\lambda_1')
                
            end
        end
        % Add legend
        legend();
        
        % Include suptitle
        suptitle(plotDescrStruct.suptitle(paramArray{parID}.plotDescr))
       
        % Save figure
        figureSavingName =  figureFolderName + "/" + ...
            plotDescrStruct.plotSavingName + "_" + ...
            simulationSetting+'_'+paramName +...
            "_N" + titleN + ...
            "_T" + titleT ;
        
        % Save both EPS and JPG versions
        saveas(p,figureSavingName,'epsc');
        saveas(p,figureSavingName,'jpg')
    end
end

%% Weights

customColors = [
    0.2 0.4 0.6;  % Light blue
    0.1 0.3 0.5;  % Medium blue
    0.0 0.2 0.4   % Dark blue
];


