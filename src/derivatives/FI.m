function [RFI, m, p] = FI(x,alpha,A,y,b,c,Q,step,d,w)
% FI compute the Fisher Information for x & alpha and calculate the
% cholesky factorization of the FI and mean of the MALA proposal
%
% [RFI, m, p] = FI(x,alpha,A,y,b,c,Q,step,d,w)
% Input:
% x         reconstructed field
% alpha     precision paramter of Dirichlet Distribution
% A         location matrix, connecting the data to the grid cells
% y         data
% b,c       hyper parameters of alpha
% Q         precision matrix of Gausian Markov Random Field 
% step      current step size of adaptive MALA
% d         d = D-1 the dimention of the transformed D-composional data,
%           here the link function is alr
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
%
% output:
% RFI       is the cholescky factorization of Fisher information
% m         mean of the MALA proposal
% p         the approximate minimum degree permutation vector for the 
%           sparse matrix, check amd for more info
%
% FI.m 2018-07-14 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

if isempty(w), w=1; end

Ax=reshape(A*(x(:)),[(size(A*x(:),1))/d ,d]);
Az = invalr(Ax);
Az_w = bsxfun(@times, w, Az);

[~, ~, ~, H] = hessian(Az,alpha,w,y);

% H is second derivatives of loglikelihood (H-> FI should multiply with -)
X2 = -A'*H*A + Q;
alpha2 = sum(-w.^2.*psi(1,alpha*w) + sum(Az_w.^2.*psi(1,alpha*Az_w),2))+(b-1)/alpha^2;

dzsj = dalr(Az);
Xalpha = sum(repmat(alpha*bsxfun(@times,w,Az_w),[d,1]).*dzsj.*repmat(psi(1,alpha*Az_w),[d,1]),2);

FI = [X2,       A'*Xalpha;
      Xalpha'*A, alpha2]; 

p = amd(FI);
FI = FI(p,p);

[RFI, err] = chol(FI);

if err~=0
  %add a small diagonal element for numerical stability
  RFI = chol(FI + 1e-3*speye(size(FI)));
end

% dL is the first derivatives of negative likelihood
[~, dlog] = dL(x,alpha,A,y,b,c,Q,d,w);
dl = -dlog;
dl = dl(p);

old = [x(:);alpha];
old = old(p);

m = old + 0.5*step^2*(RFI\(RFI'\dl));

end