function logdetQ = logDetQ(kappa,lambda)
% LOGDETQ compute logarithm of determinant of Q, log(|Q|)
% Q is precision matrix of Gaussian Markov Random Field
%
% logdetQ = logDetQ(kappa,lambda)
% input:
% kappa     the spatial range parameter of the fields
% lambda    eigen value of UDU' of matrix Q
%           computed using the function compLambda.m
%
% output:
% logdetQ   log(|Q|)
%
% logDetQ.m 2018-07-13 Behnaz@pirzamanbin.name$
% Reference https://arxiv.org/abs/1511.06417

ksi = kappa^4 + 2*kappa^2*lambda(:) + lambda(:).^2;
logdetQ = sum(log(ksi));
end
