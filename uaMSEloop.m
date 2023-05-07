% Script that computes all the MSE 

% Infeasible
mseInfeasibleFixed(:, i, : ) = mean(errInLoopInfeasiblePsiFixed.^2, 2);
mseInfeasibleFixedRandom1(:, i, : ) = mean(errInLoopInfeasiblePsiLargeRandom1.^2, 2);
mseInfeasibleFixedRandom2(:, i, : ) = mean(errInLoopInfeasiblePsiLargeRandom2.^2, 2);
mseInfeasibleFixeSmallBias1(:, i, : ) = mean(errInLoopInfeasiblePsiLargeSmallBias1.^2, 2);
mseInfeasibleFixeSmallBias2(:, i, : ) = mean(errInLoopInfeasiblePsiLargeSmallBias2.^2, 2);
mseInfeasibleFixeLargeBias1(:, i, : ) = mean(errInLoopInfeasiblePsiLargeLargeBias1.^2, 2);
mseInfeasibleFixeLargeBias2(:, i, : ) = mean(errInLoopInfeasiblePsiLargeLargeBias2.^2, 2);


% Plug-in
msePlugInFixed(:, i, : ) = mean(errInLoopPlugInFixed.^2, 2);
msePlugInLargeRandom1(:, i, : ) = mean(errInLoopPlugInLargeRandom1.^2, 2);
msePlugInLargeRandom2(:, i, : ) = mean(errInLoopPlugInLargeRandom2.^2, 2);
msePlugInLargeSmallBias1(:, i, : ) = mean(errInLoopPlugInLargeSmallBias1.^2, 2);
msePlugInLargeSmallBias2(:, i, : ) = mean(errInLoopPlugInLargeSmallBias2.^2, 2);
msePlugInLargeLargeBias1(:, i, : ) = mean(errInLoopPlugInLargeLargeBias1.^2, 2);
msePlugInLargeLargeBias2(:, i, : ) = mean(errInLoopPlugInLargeLargeBias2.^2, 2);


% AIC
mseAIC(:, i, : ) =  mean(errInLoopAIC.^2, 2);

% MMA
mseMMA(:, i, : ) =  mean(errInLoopMMA.^2, 2);

% Individual
mseIndividual(:, i) =  mean(errInLoopIndividual.^2, 2);

% MG
mseMG(:, i, : ) =  mean(errInLoopMG.^2, 2);