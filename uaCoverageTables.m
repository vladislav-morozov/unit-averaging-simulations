
parID = 4;
tableID = 2;
if tableID == 1 % coverage
    valueTable = coverageNT;
else
    valueTable = lengthNT;
end

paramName = paramArray{parID}.saveName;

numCols = length(theta1Range);

outputTable = [];
for nID=1:2
    for tID=1:2
        % Extract coverage values
        indVals = valueTable{nID, tID}.(paramName){:, 1}';
        unrestrVals = valueTable{nID, tID}.(paramName){:, 5}';
        
        % Individual values go in odd positions, fixed-N in even
        currentRow = nan(1, numCols*2);
        currentRow( (1:numCols)*2-1) = indVals;
        currentRow((1:numCols)*2) = unrestrVals;
        
        % Create table
        outputTable = [outputTable; currentRow];
        
    end
end
% Write the table
table2latex(array2table(round(outputTable,2)))