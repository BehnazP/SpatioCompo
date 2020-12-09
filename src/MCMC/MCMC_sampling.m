function MCMC = MCMC_sampling(theta0,priors,y,w,A,B,iter)
% MCMC_SAMPLING Samples a new set of alpha,beta,x,kappa and rho using a
%               metropolis adjusted langevin algorithm (MALA) for alpha,
%               beta and x and MH-gibbs step for kappa and rho given x.
% For reference see https://arxiv.org/abs/1511.06417
%
% MCMC = MCMC_sampling(theta0,priors,data,w,A,B,iter)
%
% Input
% theta0    [alpha0 (1 x 1) ,kappa0 (1 x 1) ,rho0 (d x d) ,x0 (N x d), beta0 (p x d)]
%
% priors    [alpha.a,alpha.b,...
%           kappa.range,kappa.sigma,...
%           rho.Sigma,rho.df,...
%           field.G,field.lambda,...
%           logstep.MALA,logstep.Gibbs]
%
% data      D compositional data between (0,1)
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
% A         location matrix, connecting the data to the grid cells (N x N)
% B         covariates such as (1, elevation and etc) (N x p)
% iter      number of MCMC sampling iteration
%
% Output
% MCMC      structure
%           samples [.alpha .beta .x .kappa .rho]
%           acceptance rate related [.MALA.acc .MALA.count .MALA.logstep .Gibbs.acc .Gibbs.count .Gibbs.logstep]
%
% MCMC_sampling.m 2018-06-20 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

if isempty(w), w=1; end

[N, d] = size(theta0.x);

% making A1=[B, A] & x0=[beta;x] &
AB = [B, A];
old = [theta0.beta(:);theta0.x(:)];
p = size(B,2)/d;

sigma_beta = repmat(1./gamrnd(1, 2*ones(p,1))*1/gamrnd(1,2),[2,1]);
MCMC.l2 = zeros(p, iter); MCMC.l2(:,1) = 1;
MCMC.tau2 = zeros(iter,1); MCMC.tau2(1) = 1;

MCMC.alpha = zeros(iter,1); MCMC.alpha(1) = theta0.alpha;
MCMC.kappa = zeros(iter,1); MCMC.kappa(1) = theta0.kappa;
MCMC.rho = zeros(d*d,iter);MCMC.rho(:,1) = theta0.rho(:);

if iter>1e5
    MCMC.x = zeros(N*d,(iter/50)); MCMC.x(:,1) = theta0.x(:);
else
    MCMC.x = zeros(N*d,(iter)); MCMC.x(:,1) = theta0.x(:);
end
MCMC.beta = zeros(p*d,iter); MCMC.beta(:,1) = theta0.beta(:);

MCMC.Gibbs.logstep = zeros(iter,1); MCMC.Gibbs.logstep(1) = priors.logstep.Gibbs;
MCMC.MALA.logstep = zeros(iter,1); MCMC.MALA.logstep(1) = priors.logstep.MALA;

MCMC.Gibbs.count = zeros(iter,1);
MCMC.MALA.count = zeros(iter,1);

h = waitbar(0,'Processing...','Name','MCMC Sampling');
for i = 2:iter
    waitbar(i/iter,h)

    rho = reshape(MCMC.rho(:,i-1),[d,d]);
    [~, rhoxQ] = Q_rhoxQ(rho,MCMC.kappa(i-1),priors.field.G);
    Q0 = blkdiag(speye(2*p)./sigma_beta,rhoxQ);

    %sample alpha, beta and x
    [new, MCMC.MALA.count(i), MCMC.MALA.logstep(i)] = MALA(old,MCMC.alpha(i-1),y,w,Q0,AB,priors,d,MCMC.MALA.logstep(i-1),i);

    MCMC.alpha(i) = new(end);%
    MCMC.beta(:,i) = new(1:d*p);%
    if iter>1e5
        if (mod(i,50)==0)
             MCMC.x(:,(i/50))= new(d*p+1:(N+p)*d);
        end
    else
        MCMC.x(:,i)= new(d*p+1:(N+p)*d);
    end

    old = new(1:end-1);
    tmp = new(2*p+1:(N+p)*d);
    x = reshape(tmp,[size(MCMC.x(:,i),1)/d ,d]);

    %sample kappa and rho given x
    if (theta0.kappa ~= 0)
        [samplenew, MCMC.Gibbs.count(i),MCMC.Gibbs.logstep(i)] = Gibbs_kappa_rho(x,MCMC.kappa(i-1),rho,priors,MCMC.Gibbs.logstep(i),i);
        MCMC.kappa(i) = samplenew.kappa;
        MCMC.rho(:,i) = samplenew.rho(:);
    else
        MCMC.rho(:,i) = invwishart(x,0,priors);
        MCMC.Gibbs.logstep = 0;
    end

    %----------------------------------------------------------------------------
    %sample auxiliary variables
    temp = 1+MCMC.l2(:,i-1).^(-1);
    phi = 1./gamrnd(1, 1./temp);
    temp = 1+MCMC.tau2(i-1)^(-1);
    xi = 1/gamrnd(1,1/temp);

    %sample lambda
    tempa = (1+d)/2;
    tempb = 1./phi + sum(reshape(MCMC.beta(:,i).^2,d,p))'/(2*MCMC.tau2(i-1));
    MCMC.l2(:,i) = 1./gamrnd(tempa,1./tempb);

    %sample tau
    tempa = (1+p*d)/2;
    tempb = 1/xi + sum(sum(reshape(MCMC.beta(:,i).^2,d,p))'./MCMC.l2(:,i))/2;
    MCMC.tau2(i) = 1/gamrnd(tempa, 1/tempb);

    sigma_beta = repmat(sqrt(MCMC.tau2(i)*MCMC.l2(:,i)),[2,1]);
    %----------------------------------------------------------------------------
end

MCMC.Gibbs.acc = sum(MCMC.Gibbs.count)/iter;
MCMC.MALA.acc = sum(MCMC.MALA.count)/iter;
close(h)
end
