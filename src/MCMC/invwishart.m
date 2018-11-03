function rho = invwishart(x,kappa,priors)
% INVWISHART produces an inverse wishart random matrix of size priors.rho.Sigma
%            based on conjugate prior,
%            if rho = invwishart(df,Sigma) and x|kapa,rho = GMRF(0,rho*Q(kappa))
%            then rho|x,kappa = invwishart (N+df,Sigma+x'Qx)
%
% rho = invwishart(x,kappa,priors)
% Input:
% x         a N x d vector of reconstructed spatial field
% kappa     a spatial range parameter of a spatial field
% priors    a structure containing hyper paramteres of rho and field
%           priors.field.G
%           priots.rho.df, priors.rho.Sigma
%
% output:
% rho       a d x d variance matrix of a spatial field
%
% Gibbs_kappa_rho.m 2018-07-12 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

N = size(x,1);
Q = Q_rhoxQ([],kappa,priors.field.G);

rho = iwishrnd([],priors.rho.df+N,inv(chol(priors.rho.Sigma+x'*Q*x))');

end
