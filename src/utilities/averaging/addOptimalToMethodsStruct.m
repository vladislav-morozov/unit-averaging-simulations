function preparedMethodsArray = addOptimalToMethodsStruct(methodsArray, optimalSchemes)
% addOptimalToMethodsStruct Expands a generic optimal method into specific methods.
%
% This function replaces a generic "optimal" method within a methods array
% with specialized optimal methods based on the provided schemes. Each
% scheme defines a name, description, and associated parameters for an
% optimal method.
%
% Args:
%     methodsArray (cell array): 
%         A cell array where each element is a struct with the fields:
%             - weightFunction (function handle): The weighting function.
%             - shortName (string): A short identifier for the method.
%             - longName (string): A descriptive name for the method.
%         One element must have `shortName` equal to 'opt', representing 
%         the generic optimal method to be expanded.
%     optimalSchemes (cell array): 
%         A cell array of structs, each defining a specific optimal method 
%         with the following fields:
%             - shortName (string): A short identifier for the method.
%             - longName (string): A descriptive name for the method.
%             - unrestrictedArray (function): Function that determines
%             unrestricted units.
%
% Returns:
%     preparedMethodsArray (cell array): 
%         A cell array where each element is a struct with the fields:
%             - weightFunction (function handle): The weighting function.
%             - shortName (string): A short identifier for the method.
%             - longName (string): A descriptive name for the method.
%         For optimal methods, the struct also includes:
%             - unrestrictedArray (function): Parameters copied from the 
%               corresponding scheme in `optimalSchemes`.

    % Search for the position of the optimal method and raise an error if
    % not found.
    optimalInArray = false;
    for methodID = 1:length(methodsArray)
        if methodsArray{methodID}.shortName == "opt"
            optimalIdx = methodID;
            optimalInArray = true;
            break;
        end
    end

    % Ensure the optimal method is present in the array.
    assert(optimalInArray,...
        'The "opt" method is missing from the methods array.');

    % Extract the generic optimal method.
    optimalGeneric = methodsArray{optimalIdx};

    % Remove the original optimal method from the methods array.
    preparedMethodsArray = ...
        methodsArray(1:length(methodsArray) ~= optimalIdx);
    numNonOptimal = length(preparedMethodsArray);

    % Replace the generic optimal method with specialized ones from
    % optimalSchemes.
    numOptimalSchemes = length(optimalSchemes);
    for optSchemeID = 1:numOptimalSchemes
        % Create a new optimal method using the current scheme.
        tempOptimal = optimalSchemes{optSchemeID};
        tempOptimal.weightFunction = optimalGeneric.weightFunction;

        % Append the new method to the preparedMethodsArray.
        preparedMethodsArray{numNonOptimal + optSchemeID} = tempOptimal;
    end
end
