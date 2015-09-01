% LikelihoodRatioTest.m performs a log-likelihood ratio test for the 
%    null hypothesis model (M0) and the alternative mor complex model (M).
%    These to model have to be lumpde, hence there have to exist feasible
%    parameter values for model M such that it can egnerate the behavior of
%    model M0.
%    The function returns the p-value that the null hypothesis M0 can be
%    rejected.
%
% USAGE:
% =======
% p = LikelihoodRatioTest(logL0,k0,logL,k)
%
% INPUT:
% =====
% logL0 ... log-likelihood of M0.
% k0    ... number of parameters of M0.
% logL  ... log-likelihood of M.
% k     ... number of parameters of M.
%
% OUTPUT:
% =======
% p ... p-value that the null hypothesis can be rejected.
%
% 2012/11/18

function p = LikelihoodRatioTest(logL0,k0,logL,k)

p = 1 - chi2cdf(2*(logL-logL0),k-k0);

