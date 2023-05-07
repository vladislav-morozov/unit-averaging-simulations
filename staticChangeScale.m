function thetaNew = staticChangeScale(thetaOriginal, sigma)
    % staticChangeScale Changes variance of the coefficients without
    % redrawing them

    thetaNew = sigma*thetaOriginal;

end