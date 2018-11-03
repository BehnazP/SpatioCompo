function [L,dL] = dL(x,alpha,A,y,b,c,Q,d,w)
% DL compute both loglikelihood and first derivatives of loglikelihood 
% w.r.t x and lpha
%
%
% [L,dL] = dL(x,alpha,A,y,b,c,Q,d,w)
% Input:
% x         reconstructed field
% alpha     precision paramter of Dirichlet Distribution
% A         location matrix, connecting the data to the grid cells
% y         data
% b,c       hyper parameters of alpha
% Q         precision matrix of Gausian Markov Random Field 
% d         d = D-1 the dimention of the transformed D-composional data,
%           here the link function is alr
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
%
% output:
% L         a vector of negative loglikelihood
% dL        a vector of derivatives of negiative loglikelihood
%
% DL.m 2018-07-15 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

if isempty(w), w=1; end

x1=reshape(A*(x(:)),[(size(A*x(:),1))/d ,d]);
Az = invalr(x1);
Az_w = bsxfun(@times, w, Az);
    
logpy =  gammaln(alpha.*w) - sum(gammaln(alpha*Az_w),2) + sum(((alpha*Az_w-1).*log(y)),2);
sumlogpy = sum(logpy);

logpx =  - 0.5*(x(:))'*Q*(x(:));
logpg = (b-1)*log(alpha)-c*alpha;

L = -(sumlogpy + logpx + logpg);
    
if nargout >1
    %dLx, g is the first derivatives of loglikelihood
    [g, ~] = gradient(Az,alpha,w,y);
    dLx = A'*g - Q*(x(:));

    %dLalpha
    dlogy =  w.*psi(alpha.*w) - sum(Az_w.*psi(alpha*Az_w),2) + sum(Az_w.*log(y),2);
    sumdlogy = sum(dlogy);
    dLalpha = sumdlogy + (b-1)/alpha - c;

    % derivative of loglikelihood
    dL = [dLx;dLalpha];
    %derivatives of negiative loglikelihood
    dL = -dL;
end
end
