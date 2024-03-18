% Infeasible
errInLoopInfeasiblePsiFixed = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeRandom1 = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeRandom2 = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeSmallBias1 = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeSmallBias2 = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeLargeBias1 = zeros(numParams, numReplications, numN);
errInLoopInfeasiblePsiLargeLargeBias2 = zeros(numParams, numReplications, numN);

% Plug-in
errInLoopPlugInFixed = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeRandom1 = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeRandom2 = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeSmallBias1 = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeSmallBias2 = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeLargeBias1 = zeros(numParams, numReplications, numN);
errInLoopPlugInLargeLargeBias2 = zeros(numParams, numReplications, numN);

% AIC
errInLoopAIC = zeros(numParams, numReplications, numN);

% MMA
errInLoopMMA = zeros(numParams, numReplications, numN);

% Individual
errInLoopIndividual = zeros(numParams, numReplications);

% Mean group
errInLoopMG = zeros(numParams, numReplications, numN);