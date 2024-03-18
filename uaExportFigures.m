
% Plotting parameters
plotLineThickness= 1.1;
plotW = 1600; % width in pixels
plotH = 400; % height in pixels

%% Relative and Absolute MSE

mseSuptitles = {"relative", "absolute"};
mseStartPlotting = [2, 1];
for msePlotID = 1:2
    % Loop over parameters
    for parID = 1:numParams
        
        paramName = paramArray{parID}.saveName;
        % Loop over N and T
        if plotQuietly ~= 1
            % Plot
            p = figure('Renderer', 'painters', 'Position', [50 50 plotW plotH]);
        else
            p = figure('visible','off', 'Renderer', 'painters', 'Position', [50 50 plotW plotH]);% Do not plot
        end
        
        for tID = 1:numT
            for nID = 1:numN
                ax = subplot(numT, numN, (tID-1)*numN+nID);
                currentParTable = mseTablesN{nID, tID}.(paramName);
                numApproaches = size(currentParTable, 2);
                ax.LineStyleOrderIndex = ax.LineStyleOrderIndex; % [1]
                ax.LineStyleOrder = {'-o','--+','-*','--x','-s','--d','-v','->','-h','-^'};
                %     ax.ColorOrder = [0 0 0 ];
                hold(ax, 'on');
                for approachID = mseStartPlotting(msePlotID):numApproaches
                    %                 if ~plotMethod(approachID)
                    %                    continue
                    %                 end
                    %
                    approach = methodsForPlotting{approachID};
                    approachName = approach.longName;
                    if msePlotID==1
                        % Relative MSE
                        MSE  = currentParTable{:, approachID}./currentParTable{:, 1};
                    else
                        MSE  = currentParTable{:, approachID};
                    end
                    plot(theta1Range, MSE, ...
                        'LineStyle', approach.lineStyle, ...
                        'Marker', approach.marker, ...
                        'LineWidth', plotLineThickness, ...
                        'Color', approach.color, ...
                        'DisplayName', approachName)
                    
                    
                end
                hold(ax,'off');
                xlim([min(theta1Range), max(theta1Range)])
                if msePlotID==1
                    ylim([yMinRelMSEPlot, yMaxRelMSEPlot])
                else
                    ylim([0, 1.2*max(currentParTable.unrestr)])
                end
                
                h = yline(1);
                h.HandleVisibility='off';
                
                title("N="+num2str(valuesN(nID))+", T="+num2str(valuesT(tID)))
                xlabel('\lambda_1')
                
            end
        end
        legend();
        if msePlotID==1
           % Absolute MSE
           suptitle('Averaging estimator, ' + ...
               paramArray{parID}.plotDescr + ...
               ', ratio of MSE to individual estimator')
        else
            % Relative MSE
            suptitle('Averaging estimator, ' + ...
                paramArray{parID}.plotDescr + ...
                ', MSE')
        end
         
        % Save figure
        figureSavingName =  figureFolderName + "/" + ...
            mseSuptitles{msePlotID} + "_" + ...
            simulationSetting+'_'+paramName +...
            "_N" + titleN + ...
            "_T" + titleT ;
        
        % Save both EPS and JPG versions
        saveas(p,figureSavingName,'epsc');
        saveas(p,figureSavingName,'jpg')
    end
end

%% Something
