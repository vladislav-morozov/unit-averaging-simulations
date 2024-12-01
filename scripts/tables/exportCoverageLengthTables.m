% ===========================================================
% File: exportCoverageLengthTables.m
% Description: This script exports the figures based on simulations for the
%              mean squared error, bias, variance, weights, and
%              unrestricted units for unit averaging.
%
% ===========================================================
%
% Project Name: Unit Averaging for Heterogeneous Panels
% Authors: Christian Brownlees, Vladislav Morozov
%
% Created tables are saved in latex formate in the 'results/tables' folder.
% ===========================================================

% Load current file 
fileSaveName = makeOutputFileName('CI', simulationSetting, ...
                                  numReplicationsCI, valuesN, valuesT);
load(fileSaveName)

% Determine number of columns: each column corresponds to a grid point
numCols = length(theta1Range);

% Compute number of sample sizes used
numN = length(valuesN);
numT = length(valuesT);

% Loop through parameters
for parID = 1 : length(paramArray)

    % Extract parameter name
    paramName = paramArray{parID}.saveName;

    % Loop through table types (coverage and length)
    for tableID = 1 : 2
        % Determine table type 
        if tableID == 1 
            % Coverage
            valueTable = coverageNT;       % data source
            propertyName = 'coverage';     % part of table file name
        else
            % Length
            valueTable = lengthNT;    
            propertyName  = 'length'; 
        end

        % Create array to hold the table
        outputTable = [];

        % Create array to hold row and column labels
        rowLabels = cell(numN * numT, 1);
        colLabels = cell(numCols * 2, 1);

        % Loop through the  values of (N, T)
        for nID = 1 : numN 
            for tID = 1 : numT

                % Extract coverage/length values
                indVals = valueTable{nID, tID}.(paramName){:, 1}';
                unrestrVals = valueTable{nID, tID}.(paramName){:, 5}';

                % Individual values go in odd positions, fixed-N in even
                currentRow = nan(1, numCols * 2);
                currentRow((1 : numCols) * 2 - 1) = indVals;
                currentRow((1 : numCols) * 2) = unrestrVals;

                % Expand table with current row
                outputTable = [outputTable; currentRow];
                
                % Create the current row label
                rowLabels{(nID - 1)*numT + tID} = ...
                    sprintf('N=%d, T=%d', valuesN(nID), valuesT(tID));
            end
        end
        % Create column labels
        colLabels((1 : numCols) * 2 - 1) = ...
            arrayfun(@(gridPoint) ...
            sprintf("Grid point %0.2f, Individual", gridPoint), ...
            theta1Range, 'UniformOutput',false);
        colLabels((1 : numCols) * 2 ) = ...
            arrayfun(@(gridPoint) ...
            sprintf("Grid point %0.2f, Fixed-N", gridPoint), ...
            theta1Range, 'UniformOutput',false);

        % Recast the data array into a table 
        outputTable = array2table(round(outputTable,2), ...
            "RowNames",rowLabels, ...
            'VariableNames',[colLabels{:}]);

        % Write the table
        tableSavingName = makeTableName(...
            simulationSetting, ...
            propertyName, paramName, ...
            valuesN, valuesT);
        table2latex(outputTable, tableSavingName)
    end
end 

%% Auxiliary function: table name generator
function tableSavingName = makeTableName(...
    simulationSetting, ...
    propertyName, paramName, ...
    valuesN, valuesT)
    % makeTableName Generates a standardized table name for export
    %
    % Args:
    %     leadingString (string): String to insert at the beginning of the
    %         file name.
    %     simulationSetting (string): Name of the data-generating process.
    %     propertyName (string): Coverage or length.
    %     paramName (string): Name of current target parameter.
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
    tableSavingName = sprintf('results/tables/%s_%s_%s_N%s_T%s.tex', ... 
        simulationSetting, ...      % Simulation context
        propertyName, ...           % Coverage or length
        paramName, ...              % Parameter name
        titleN, titleT);            % Sample sizes and cross-sections 
end