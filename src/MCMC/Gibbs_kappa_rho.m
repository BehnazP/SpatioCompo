function [sample, count, logstep2] = Gibbs_kappa_rho(x,kappa,rho,priors,logstep0,iter)
% GIBBS_KAPPA_RHO an adaptive Metropolis adjusted (Reimann manifold) Langevin is used 
%                 to propose a new value for x and alpha
%
% [sample, count, logstep2] = Gibbs_kappa_rho(x,kappa,rho,priors,logstep0,iter)
% Input:
% x         a N x d vector of reconstructed spatial field
% kappa     a spatial range parameter of a spatial field
% rho       a d x d variance matrix of a spatial field
% priors    a structure containing hyper paramteres of kappa and rho
% logstep0  logarithm of the current variance of the adaptive random walk algorithm
% iter      current iteration in MCMC sampling
%
% output:
% sample    structure contain:
%           .kappa proposed range
%           .rho proposed variances
% count     indication of the acceptance of MH proposal;
%           1 if the proposed MH is accpeted and
%           0 if the old value remains
% logstep2  new value for the logarithm of the varaince of the adaptive random walk
%
% Gibbs_kappa_rho.m 2018-07-12 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

step = exp(logstep0/2);

if kappa <= 0
	acc_prob = 0;
else
    old = kappa;
    prop = log(old) + step*randn(1);
    new = exp(prop);
    logpxy = -log_posterior_kappax(x,old,priors);
    logpyx = -log_posterior_kappax(x,new,priors);
    acc_prob = exp(min(0,logpyx-logpxy+log(new)-log(old)));
end

U = rand(1);
if (U < acc_prob)
    sample.kappa = new ;
    sample.rho = invwishart(x,new,priors);
    count = 1;
else
    sample.kappa = kappa;
    sample.rho = rho;
    count = 0;
end

logstep2 = logstep0 + (iter+1)^(-0.5)*(acc_prob - 0.44);

end
