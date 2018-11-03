function [point_dist, mean_dist] = compo_dist(data1,data2)
% compo_dist compute the distance between two D-compositional data 
%
% [point_dist, mean_dist] = compo_dist(data1,data2)
% 
% data1 and data2 are two D-compositional data
%
% dist(data1,data2) = ([alr(data1)-alr(data2)]^t*H*[alr(data1)-alr(data2)])^1/2
% where H = [2, 1;
%            1, 2];
%
% Return point_dist which is point wise distace between two vectors and 
% mean_dist which is the mean distance
% 
% compo_dist.m 2018-07-13 Behnaz@pirzamanbin.name$
% Reference 
% B. Pirzamanbein (2018) https://arxiv.org/abs/1511.06417
% J. Aitchison (2000) https://link.springer.com/article/10.1023%2FA%3A1007529726302


if size(data1,1)~=size(data2,1)
   error('Bad size for data')
end

f = data1./data2;
g =  bsxfun(@rdivide,f,sum(f,2));

D = size(data1,2);
H = eye(D-1)+ones(D-1);

dist = diag((alr(g))*(H\(alr(g))'));
point_dist = sqrt(dist);
mean_dist = mean(point_dist);

end