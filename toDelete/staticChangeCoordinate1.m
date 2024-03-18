function betaNew = staticChangeCoordinate1(b1, beta)
% staticChangeCoordinate1 Replaces first coordinate of coefficient vector
% beta with new value b1

    betaNew = beta;
    betaNew(1, 1) = b1;
end