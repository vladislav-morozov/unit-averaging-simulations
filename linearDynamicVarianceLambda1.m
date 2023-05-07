function var1 = linearDynamicVarianceLambda1(lambda, beta,sigmaSq)

    var1 = sigmaSq*((1-(lambda.^2))./((beta.^2+sigmaSq)));
end