function [sample, count, logstep2] = MALA(x,alpha,y,w,Q,A,priors,d,logstep0,iter)
% MALA an adaptive Metropolis adjusted (Reimann manifold) Langevin is used
% to propose a new value for x and alpha
%
% [sample, count, logstep2] = MALA(x,alpha,y,w,Q,A,priors,d,logstep0,iter)
% Input:
% x         reconstructed field
% alpha     precision paramter of Dirichlet Distribution
% y         data
% w         location related weigths if they exist, i.e. if there is some
%           probability based on the location of the data
% Q         precision matrix of Gausian Markov Random Field
% A         location matrix, connecting the data to the grid cells
% priors    structure containing:
%           alpha.a,alpha.b hyper parameters of alpha
% logstep0  logarithm of current step size of adaptive MALA
% iter      current iteration in MCMC sampling
%
% output:
% sample    proposed x and alpha
% count     indication of the acceptance of MALA proposal;
%           1 if the proposed MALA is accpeted and
%           0 if the old value remains
% logstep2  new value for logarithm of step size of MALA
%
% MALA.m 2018-07-12 writen by Johan Lindstrom and Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

step = exp(logstep0/2);
old = [x(:);alpha];

[RFI0, m0, p0] = FI(old(1:end-1),old(end),A,y,priors.alpha.a,priors.alpha.b,Q,step,d,w);
epsilon = randn(size(RFI0,1),1);
new = m0 + RFI0\(step*epsilon);
new(p0) = new;
alphanew = new(end);

if (alphanew<0)
    acc_prob = 0;
else
    logpxy = dL(old(1:end-1),old(end),A,y,priors.alpha.a,priors.alpha.b,Q,d,w);
    logpyx = dL(new(1:end-1),new(end),A,y,priors.alpha.a,priors.alpha.b,Q,d,w);

    if isfinite(logpyx)

        [RFI, m, p] = FI(new(1:end-1),new(end),A,y,priors.alpha.a,priors.alpha.b,Q,step,d,w);

        old = old(p);
        %and proposal density: q(old|new)
        %the above is simplified since:
        %   (old-m)' * (FI/step^2) * (old-m) =
        % = (old-m)' * (RFI/step)' * RFI/step * (old-m) =
        % = ||RFI/step * (old-m)||_2^2
        logqxy = sum(log(diag(RFI))) - sum(((RFI./step)*(old-m)).^2)/2;

        %and proposal density: q(new|old)
        %the above is simplified since:
        %   (new-m)' * (FI_0/step^2) * (new-m) =
        % = (RFI0 \ step*epsilon)' (FI0/step^2) * (RFI0 \ step*epsilon) =
        % = epsilon' * inv(RFI0)' * FI0 *inv(RFI0) * epsilon) =
        % = epsilon' * epsilon
        logqyx = sum(log(diag(RFI0))) - sum(epsilon.^2)/2;

        acc_prob = exp( min(0,-logpyx+logpxy+logqxy-logqyx));
    else
        acc_prob = 0;
    end
end

U = rand(1);
if (U < acc_prob)
    sample = new;
    count = 1;
else
    sample = [x(:);alpha];
    count = 0;
end

logstep2 = logstep0 + (iter+1)^(-0.5)*(acc_prob - 0.56);

end
