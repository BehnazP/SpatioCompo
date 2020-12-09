function transform = invalr(compo)
% INVALR is inverse alr transformation, i.e. inverse additive log ratio
%
% transform = invalr(compo)
% input:
% compo         a n x d data between interval (-Inf +Inf)
%
% Output:
% transform     a n x D vector of compositional data (D = d+1) with properties;
%               1- sum(compo) = 1, and
%               2- each component of compo is belong to the interval (0,1)
%
% INVALR.m 2018-07-13 Behnaz@pirzamanbin.name$

[n, d] = size(compo);
D = d+1;
transform = zeros(n,D);

for k = 1:d
    transform(:,k)=exp(compo(:,k))./(1 + sum(exp(compo),2));
end
transform(:,D)=1./(1 + sum(exp(compo),2));
end
