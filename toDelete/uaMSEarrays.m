% Infeasible
mseInfeasibleFixed = zeros(numParams, etaRangeLen, numN);
mseInLargeRandom1 = zeros(numParams, etaRangeLen, numN);
mseInLargeRandom2 = zeros(numParams, etaRangeLen, numN);
mseInfeasibleLargeSmallBias1 = zeros(numParams, etaRangeLen, numN);
mseInfeasibleLargeSmallBias2 = zeros(numParams, etaRangeLen, numN);
mseInfeasibleLargeLargeBias1 = zeros(numParams, etaRangeLen, numN);
mseInfeasibleLargeLargeBias2 = zeros(numParams, etaRangeLen, numN);

% Plug-in
msePlugInFixed = zeros(numParams, etaRangeLen, numN);
msePlugInLargeRandom1 = zeros(numParams, etaRangeLen, numN);
msePlugInLargeRandom2 = zeros(numParams, etaRangeLen, numN);
msePlugInLargeSmallBias1 = zeros(numParams, etaRangeLen, numN);
msePlugInLargeSmallBias2 = zeros(numParams, etaRangeLen, numN);
msePlugInLargeLargeBias1 = zeros(numParams, etaRangeLen, numN);
msePlugInLargeLargeBias2 = zeros(numParams, etaRangeLen, numN);

% AIC
mseAIC = zeros(numParams, etaRangeLen, numN);

% MMA
mseMMA = zeros(numParams, etaRangeLen, numN);

% Individual
mseIndividual = zeros(numParams, etaRangeLen);

% Mean group
mseMG = zeros(numParams, etaRangeLen, numN);