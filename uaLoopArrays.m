% Infeasible
errInLoopInfeasiblePsiFixed = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeRandom1 = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeRandom2 = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeSmallBias1 = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeSmallBias2 = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeLargeBias1 = zeros(numPar, numReplications, Nlen);
errInLoopInfeasiblePsiLargeLargeBias2 = zeros(numPar, numReplications, Nlen);

% Plug-in
errInLoopPlugInFixed = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeRandom1 = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeRandom2 = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeSmallBias1 = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeSmallBias2 = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeLargeBias1 = zeros(numPar, numReplications, Nlen);
errInLoopPlugInLargeLargeBias2 = zeros(numPar, numReplications, Nlen);

% AIC
errInLoopAIC = zeros(numPar, numReplications, Nlen);

% MMA
errInLoopMMA = zeros(numPar, numReplications, Nlen);

% Individual
errInLoopIndividual = zeros(numPar, numReplications);

% Mean group
errInLoopMG = zeros(numPar, numReplications, Nlen);