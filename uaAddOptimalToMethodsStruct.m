function preparedMethodsArray = ...
    uaAddOptimalToMethodsStruct(methodsArray, optimalSchemes)
% uaAddOptimalToMethodsStruct Expands the generic Optimal method in
% methodsArray into specific optimal methods with names and unrestricted
% patters determined by the entries of the optimalSchemes array
%
% Inputs: 1. methodsArray -- cell array. Each element is a struct
%          with fields .weightFunction, .shortName, .longName. One element 
%          must have shortName 'opt'. This element will be expanded using 
%          elements of optimalSchemes
%         2. optimalSchemes -- cell array. Each element is a struct with
%         fields .shortName, .longName, and .unrestrictedArray. Each 
%         element will be inserted as an optimal approach into 
%         preparedMethodsArray
%
% Returns: preparedMethodsArray -- cell array. Each element is a struct
%          with fields .weightFunction, .shotName, and .longName. Optimal 
%          methods also have third field .unrestrictedArray, taken from the
%          corresponding field in optimalSchemes.
%          The .weightFunction has the signature of the .weightFunctions in
%          methodsArray. The names of the new optimal methods is determined
%          by the optimalSchemes names. 

% Search for the position of optimal method and raise error if not found
optimalInArray = false;
for methodID = 1:length(methodsArray)
    if methodsArray{methodID}.shortName == "opt"
       optimalIdx = methodID;
       optimalInArray  = true;
       break 
    end
end

% Raise error if Optimal weights not in the array
assert(optimalInArray, 'Optimal averaging not in the methods array')

% Extract the generic optimal method
optimalGeneric = methodsArray{optimalIdx};

% Drop the original optimal method
preparedMethodsArray = methodsArray(1:length(methodsArray)~=optimalIdx);
numNonOptimal  = length(preparedMethodsArray);

% Replace the generic optimal method with specialized ones that carries
% the corresponding restricted units
numOptimalSchemes = length(optimalSchemes);
for optSchemeID = 1:numOptimalSchemes
    % Copy the description of the new scheme into a temp struct
    tempOptimal = optimalSchemes{optSchemeID};
    % Add the generic weight function
    tempOptimal.weightFunction = optimalGeneric.weightFunction;
    % Insert into the prepared methods array
    preparedMethodsArray{numNonOptimal+optSchemeID} = ...
        tempOptimal;
end

end