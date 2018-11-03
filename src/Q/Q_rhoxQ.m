function [Q, rhoxQ] =  Q_rhoxQ(rho,kappa,G)
% Q_RHOXQ compute a precison matrix of Gaussian Markov Random field and a
% kronecker product of a variance matrix and Q
% 
% [Q, rhoxQ] =  Q_rhoxQ(Rho,kappa,G)
% input:
% kappa     a spatial range parameter of a spatial field
% rho       a variance matrix of a spatial field (d x d)
% G         a precision matrix for laplace^(alpha/2) u = e on a regular 
%           grid with Neuman boundaries. Computed using function createQ.m
%
% output:
% Q         a precision matrix of Gaussian Markov Random field
% rhoxQ     a Nd x Nd matrix; this is a kronecker product of 
%           rho (d x d) and Q (N x N) 
% Example:  if d = 2 then
%           [rho(1,1)*Q      rho(1,2)*Q
%           rho(2,1)*Q      rho(2,2)*Q]
%
% logDetQ.m 2018-07-13 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

Q = kappa^4*speye(size(G)) + 2*kappa^2*G + G'*G;
    
if ~isempty(rho)
    rhoxQ = kron(inv(rho),Q);
end

end