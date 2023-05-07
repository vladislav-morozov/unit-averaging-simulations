function betaNew = linearDynamicDrawBimodal(a, b, sigma1, sigma2, p, N, seedCoef)
%   linearDynamicDrawBimodal Draw a new set of betas from a mixture of
%   normal N(a, sigma1) and N(b, sigma2) distributions with proportion p
%   going towards N(a, sigma1)

   rng(seedCoef,'philox')
   r = rand(N,1)>=p;
   R1 = normrnd(a,sigma1,N,1);
   R2 = normrnd(b,sigma2,N,1);
   betaNew = (1-r).*R1+r.*R2;
   betaNew(1) = a;
end