function [g, dg]= gradient(z,alpha,w,y)
% GRADIENT of log Dirichle distribution with alr 
% transformation of D-compositional data
%
% [g, dg]= gradient(z,alpha,w,y)
% input: 
% z         a n x D compositional data with properties;
%            1- sum(compo) = 1, and
%            2- each component of compo is belong to the interval (0,1)
% alpha     precision parameter of Dirichlet distribution
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
% y         data
%
% Gradient of Dirichlet (alpha * z) with respect to z
% dg = dDir/dz              n x D dimension
% g  = sum(dDir/dz*dz/dx)   n x 1  
% dz/dx = dalr(z); or any other derivative of a transformation e.g. clr
%
% output:
% g         a vector of gradient of dirichlet distribution w.r.t
%           transformed D-compositional field, here is alr trasnformed
%           
% dg        derivatives of Dirichlet distribution w.r.t D-compositional
%           field z
%
% GRADIENT.m 2018-07-14 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417 

D = size(z,2);
d = D - 1;
n = size(y,1);

dz = dalr(z);
z_w = bsxfun(@times, w, z);

dg=zeros(n,D);

for l = 1:D
    dg(:,l) = -alpha*w.*psi(alpha*z_w(:,l)) + alpha*w.*log(y(:,l));
end
g = sum(repmat(dg,[d,1]).*dz,2);    
end

