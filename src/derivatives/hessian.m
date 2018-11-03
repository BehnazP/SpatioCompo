function [d2gii, d2gij,  H, H_zgz] = hessian(z,alpha,w,y)
% HESSIAN of log Dirichlet distribution with alr transformation of
% D-composional data 
%
% [d2gii, d2gij,  H, H_zgz] = hessian(z,alpha,w,y)
% input: 
% z         a n x D compositional data with properties;
%            1- sum(compo) = 1, and
%            2- each component of compo is belong to the interval (0,1)
% alpha     precision parameter of Dirichlet distribution
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
% y         data
%
% Hessian or second derivaties of Dirichlet (alpha * z) with respect to z
% d2g = d2Dir/dz2              
% H  = d2Dir/dz2*dz/dx          
% dz/dx = dalr(z); or any other derivative of a transformation e.g. clr
%
% output: 
% d2gii   second derivatives of the Dirichlet distribtion w.r.t same component  
% d2gij   second derivatives of the Dirichlet distribtion w.r.t other components 
% H       Hessian of the Dirichlet distribution
% H_zgz   Hessian of the Dirichlet distribution without adding other
%         components
%
% HESSIAN.m 2018-07-14 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417 


[n, D] = size(z);
d = D - 1;

dz = dalr(z);
z_w = bsxfun(@times, w, z);

d2g = zeros(n,D);

for k = 1:D
    d2g(:,k) = - alpha^2*w.^2.*psi(1,alpha*z_w(:,k));
end

bigd2g = repmat(d2g,[d,1]);
d1 = bigd2g.*dz.^2;
d2 = sum (d1,2);

[~ , dg]= gradient(z,alpha,w,y);
Adg = dg;
bigdg = repmat(Adg,[d,1]);

for k = 1:d
    Z1((k-1)*n+1:(k*n),1) =  ones(n,1)-2*z(:,k);
end

C = repmat(Z1,[1,D]);
dzC = dz.*C;
d3 = sum(bigdg.*dzC,2);

d4 = d2+d3;
d2gii = reshape(d4,[n,d]);
d2gii_zgz = reshape(d2,[n,d]);

h0 = dz(1:n,:).*dz(n+1:2*n,:);
h1 = repmat(h0,[d,1]);

d5 = sum(h1.*bigd2g,2);

Z0 = [-z(:,1).*z(:,2).*Z1(1:n,:) -z(:,1).*z(:,2).*Z1(n+1:2*n,:) 2*z(:,1).*z(:,2).*z(:,3)];
Z2 = repmat(Z0,[d,1]);

d6 = sum(Z2.*bigdg,2);
d7 = d5+d6;
d2gij = reshape(d7,[n,d]);
d2gij_zgz = reshape(d5,[n,d]);

G =[d2gij(:,1)  d2gii(:,1) zeros(n,1);
    zeros(n,1)  d2gii(:,2) d2gij(:,1)];
H = spdiags(G,[-n 0 n ],n*d,n*d);

G_zgz =[d2gij_zgz(:,1) d2gii_zgz(:,1) zeros(n,1);
        zeros(n,1)     d2gii_zgz(:,2) d2gij_zgz(:,1)];
H_zgz = spdiags(G_zgz,[-n 0 n ],n*d,n*d);

end

