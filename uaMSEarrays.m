% Infeasible
mseInfeasibleFixed = zeros(numPar, endOfRange, Nlen);
mseInLargeRandom1 = zeros(numPar, endOfRange, Nlen);
mseInLargeRandom2 = zeros(numPar, endOfRange, Nlen);
mseInfeasibleLargeSmallBias1 = zeros(numPar, endOfRange, Nlen);
mseInfeasibleLargeSmallBias2 = zeros(numPar, endOfRange, Nlen);
mseInfeasibleLargeLargeBias1 = zeros(numPar, endOfRange, Nlen);
mseInfeasibleLargeLargeBias2 = zeros(numPar, endOfRange, Nlen);

% Plug-in
msePlugInFixed = zeros(numPar, endOfRange, Nlen);
msePlugInLargeRandom1 = zeros(numPar, endOfRange, Nlen);
msePlugInLargeRandom2 = zeros(numPar, endOfRange, Nlen);
msePlugInLargeSmallBias1 = zeros(numPar, endOfRange, Nlen);
msePlugInLargeSmallBias2 = zeros(numPar, endOfRange, Nlen);
msePlugInLargeLargeBias1 = zeros(numPar, endOfRange, Nlen);
msePlugInLargeLargeBias2 = zeros(numPar, endOfRange, Nlen);

% AIC
mseAIC = zeros(numPar, endOfRange, Nlen);

% MMA
mseMMA = zeros(numPar, endOfRange, Nlen);

% Individual
mseIndividual = zeros(numPar, endOfRange);

% Mean group
mseMG = zeros(numPar, endOfRange, Nlen);