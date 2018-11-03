function [Q, W, cQ] = createQ(sz, alpha, a, indW, Neumann)
% CREATEQ Creates a Q matrix for grided data of size sz.
%
%  [Q, W, cQ] = createQ(sz, alpha=2, a=1, indW=false, Neumann=true)
%
%  sz - a vector containing the size of the grid, in 1D, 2D, or higher,
%       eg. [d1 d2 d3 ...]
%  alpha - order of the laplace operator on the grid (1 or 2)
%  a - vector of same length as sz with positive anisotropy weights for
%      each increment.
%  indW - true/false return a cell-array of W for increments in each
%         direction or a single W as W = sum W_i.
%  Neumann - true/false use Neumann or torus boundary conditions.
%
%  Returns a precision matrix, Q, for
%    laplace^(alpha/2) u = e
%  on a regular grid with Neuman/torus boundaries; also returns the increment
%  matrix, W, and the circulant vector, cQ, for the matrix with torus
%  boundary conditions (can be used as a preconditioner).
%  The computed matrices assume a grid spacing of h=1, scaling for correct
%  Q is given by Q = Q*h^(alpha*(d-2))
%
%  The a-vector allows for axis anisotropy giving the Q matrix as
%    Q = ( sum( a(i) * W{i} ) )^(alpha/2)
%
% $Id: createQ.m 3724 2014-07-17 16:07:09Z johanl $

if nargin<2 || isempty(alpha), alpha=2; end
if nargin<3 || isempty(a), a=1; end
if nargin<4 || isempty(indW), indW=false; end
if nargin<5 || isempty(Neumann), Neumann=true; end

if alpha~=1 && alpha~=2
  error('alpha should be 1 or 2!')
end

if ~isscalar(a) && length(a)~=length(sz)
  error('a must be scalar or have length(sz) elements')
end
if any(a<0), error('a must be non-negative'); end
%repeat a to match length of sz.
if isscalar(a), a=repmat(a, [1 length(sz)]); end

%remove singelton dimensions
a = a(sz~=1);
sz = sz(sz~=1);

%create W matrices for each dimension
W = cell(1, length(sz));
for i=1:length(sz)
  W{i} = createW(sz(i), Neumann);
end
%Use kronecker to expand each increment along a the other dimension(s)
if length(sz)>1
  W{1} = kron( speye(prod(sz(2:end))), W{1} );
  for i=2:(length(sz)-1)
    W{i} = kron( speye(prod(sz((i+1):end))), ...
                 kron(W{i}, speye(prod(sz(1:(i-1))))) );
  end
  W{end} = kron( W{end}, speye(prod(sz(1:end-1))) );
end

%accumulate into a Q matrix
Q = a(1)*W{1};
for i=2:length(sz)
  Q = Q+a(i)*W{i};
end

%return a single W matrix rather than the above cell structure
if ~indW && nargout>1,  W = Q; end

%should Q be W (alpha=1) or W*W (alpha=2)
if alpha==2, Q = Q*Q; end

%compute circulant approximation (used as preconditioner)
if nargout>2
  d = length(sz);
  if size(sz,1)~=1, sz=sz'; end
  sz = [1 sz];
  I = [1 cumprod(sz(1:end-1))+1 cumprod(sz(1:end-1)).*(sz(2:end)-1)+1];
  cQ = sparse([2*d -ones(1,2*d)], 1, I, 1, prod(sz));

  %circulant (convolution if alpha=2)
  if alpha==2
    cQ = real(ifft( fft(full(cQ)).^2 ));
    cQ(abs(cQ) < 10*eps) = 0;
  end
end

function W = createW(n, Neumann)
W = spdiags(repmat([-1 2 -1],[n 1]),-1:1, n, n);
if Neumann
  W(1, 1) = 1; W(n, n) = 1;
else
  W(1, n) = -1; W(n, 1) = -1;
end
