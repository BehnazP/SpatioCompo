% MCMC
if isempty(iter)
    iter = 1e5;
end

priors.rho.df = 10;
priors.rho.Sigma = speye(d);

alpha0 = 10;%a start value for alpha
%hyperparameter for alpha
priors.alpha.a = 1.5;
priors.alpha.b = 0.1;
%hyper parameter for beta
sigma_beta = 1e3;
%start value for kappa and rho
kappa_start = 0;
rho_start = speye(d);

% PC prior for hyper parameter for kappa
priors.kappa.range = 1;
priors.kappa.sigma= log(100)*1/sqrt(8); % range0 = 1 there will be no range smaller than 1 implies kappa less than 1.6

priors.field.G = createQ(sz,1);
priors.field.lambda = compLambda(sz);

if ~exist('w','var')
   w=[];
end
priors.logstep.MALA = 0; % log(step^2)
priors.logstep.Gibbs = 0;

%start values for Beta and X and given Beta(1) and X(1) finding a proper
%satrt value for alpha
[~, rhoxQ] = Q_rhoxQ(rho_start,kappa_start,priors.field.G);
Q = blkdiag(speye(2*p)/sigma_beta,rhoxQ);

mu0 = [zeros(p,d);zeros(prod(sz),d)];
opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton','GradObj','on', ...
                'Hessian', 'off', 'Display', 'off','TolX',1e-1,'TolFun',1e-1);
A = [B_Rev,A_x];
betax_mode = fminunc(@(X) dL(X(1:end-1),X(end),A,y_Rev, priors.alpha.a, priors.alpha.b, Q,d,w), [mu0(:);alpha0], opts);

beta_start = reshape(betax_mode(1:p*d),[p,d]);
x_start = reshape(betax_mode(p*d+1:end-1),[prod(sz),d]);
alpha_start = betax_mode(end);

% given the mode find the proper kappa and rho start value
kappa_start = fminbnd(@(kappa) log_posterior_kappax(x_start,kappa_start,priors), 0, 2);
rho_start = invwishart(x_start,kappa_start,priors);

theta0.alpha = alpha_start;
theta0.kappa = kappa_start;
theta0.rho = rho_start;
theta0.x = x_start;
theta0.beta = beta_start;

tic;
MCMC = MCMC_sampling(theta0,priors,y_Rev,w,A_x,B_Rev,iter);
MCMC.t = toc;

%% flag_save the MCMC results

if ~isempty(flag_save)
    filename = ['MCMC_',model,'_',num2str(time),'.mat'];
    save(fullfile(pathname,filename),'MCMC')
end
