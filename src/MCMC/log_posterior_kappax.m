function logp = log_posterior_kappax(x,kappa,priors)
% LOG_POSTERIOR_KAPPAX computes the negative log likelihood of kappa given x
% 
% logp = log_posterior_kappax(x,kappa,priors)
% 
% Input:
% x         a sampled spatial field (N x d)
% kappa     a spatial range parameter of a spatial field
% priors    structure containing:
%           .kappa.sigma .kappa.range hyper parameter of kappa 
%           .rho.Sigma, .rho.df hyper parameter of rho
%           .field.G, .field.lambda parameters of field
%
% The joint posterior of kappa and rho given x is 
% p(kappa,rho|x) = p(rho|kappa,x) . p(kappa|x)
% due to conjugate prior for rho the conditionap posterior of rho is
% inverse wishahrt
% p(rho|kappa,x) = IW(Sigma+xtQx,N+df)
% therefor one can integrate out rho from the joint posterior and 
% p(kappa|x) = int( P(rho|kappa,x)P(kappa) drho ) 
%            = |Q|^d/2*P(kappa) / (|Sigma+xtQx|^(n+df)/2)
%
% output
% logp      negative loglikelihood of P(kappa|x)
%
% log_posterior_kappax.m 2018-06-20 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

[N, d] = size(x);

Q = Q_rhoxQ([],kappa,priors.field.G);
logdetQ = logDetQ(kappa,priors.field.lambda);

IxtQx = (priors.rho.Sigma + ((x)'*Q*(x)));%/(N+df);
% this give a dxd matrix in case of d = 3 the determinant is easy to calculate
logdet = 2*sum(log(diag(chol(IxtQx)))); %log(det(IxtQx));

logkappa = log(kappa);
pkappa = (priors.kappa.range-1)*logkappa - priors.kappa.sigma*kappa;

logp = -(0.5*d*logdetQ - 0.5*(N+priors.rho.df)*logdet + pkappa);
end