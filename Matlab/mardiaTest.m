function [b1p,b2p,zSkew,pSkew,zKurt,pKurt] = mardiaTest(X)
%-------------------------------------------------------------
% Mardia's Multivariate Normality Test
%
% INPUT
%   X : n × p matrix
%
% OUTPUT
%   b1p   : multivariate skewness
%   b2p   : multivariate kurtosis
%   zSkew : skewness statistic
%   pSkew : skewness p-value
%   zKurt : kurtosis statistic
%   pKurt : kurtosis p-value
%
% Reference:
% Mardia (1970)
%-------------------------------------------------------------

[n,p] = size(X);

mu = mean(X);

S = cov(X);

Sinv = pinv(S);

Xc = X-mu;

D = Xc*Sinv*Xc';

%% Mardia skewness

b1p = sum(D(:).^3)/(n^2);

chiSkew = n*b1p/6;

df = p*(p+1)*(p+2)/6;

pSkew = 1-chi2cdf(chiSkew,df);

zSkew = chiSkew;

%% Mardia kurtosis

di = diag(D);

b2p = mean(di.^2);

Expected = p*(p+2);

Var = 8*p*(p+2)/n;

zKurt = (b2p-Expected)/sqrt(Var);

pKurt = 2*(1-normcdf(abs(zKurt)));

end