function dhkj = dalr(z)
% DALR is derivative of additive log ratio transformation 
%
% dhkj = dalr(z) 
% input: 
% z         a n x D compositional data with properties;
%           1- sum(compo) = 1, and
%           2- each component of compo is belong to the interval (0,1)
%
% z is D compositinal data and
% zk = exp(xk)/1+sum(exp(xk))  if k = 1:D-1
%            1/1+sum(exp(xk))  if k = 1:D  
% then the derivatives are computed as
% dhkj = dzk/dxj = zk(1-zk) if k = j
%                  -zkzj    if k ~= j
%
% output:
% dhkj      a matrx of derivaties
%
% ALR.m 2018-07-13 Behnaz@pirzamanbin.name$

[n, D] = size(z);
d = D - 1;
dhkj=zeros(n*d,D);

    for k = 1:d
        for j = 1:D
            dhkj(1+(k-1)*n:k*n,j) = z(:,k).*(not(k-j)*ones(n,1)-z(:,j));
        end 
    end

end