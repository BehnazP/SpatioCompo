function transform = alr(compo)
% ALR is additive log ratio transformation 
%  
% transform = alr(compo)
% input: 
% compo         a n x D compositional data with properties;
%               1- sum(compo) = 1, and
%               2- each component of compo is belong to the interval (0,1)
% output:
% transform     a n x d vector of transpform composition (d = D-1) 
%               belong to interval (-Inf +Inf)
%
% ALR.m 2018-07-13 Behnaz@pirzamanbin.name$

[n, D] = size(compo);
d = D-1;
transform = zeros(n,d);
for k = 1:d
    transform(:,k)=log(compo(:,k)./compo(:,D));
end 

end