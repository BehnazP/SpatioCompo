function Lambda = compLambda(sz, a, Neumann)
% COMPLAMBDA Computes the eigenvalues of a Laplace operator on a regular grid.
%
%  Lambda = compLambda(sz, a=1, Neumann=true)
%
%  sz - a vector containing the size of the grid, in 1D, 2D, or higher, 
%       eg. [d1 d2 d3 ...]
%
%  Returns a matrix of size sz containing the eigenvalues of
%    sum( a(i) * W{i} )
%  where W{i} is the finite difference [-1 2 -1] along the i:th dimensions.
%  The eigenvalues can be used to compute W*x with DCT/FFT
%
%  1D - Example: 
%   x = randn(13,1);
%   W = createQ(length(x),1);
%   L = compLambda(length(x));
%   plot([idct(L.*dct(x)) W*x])
%
%  3D anisotropic - Example: 
%   sz = [13 7 11];
%   x = randn(sz);
%   a = rand(length(sz),1)+1;
% 
%   W_dct = createQ(sz, 1, a, [], true);
%   L_dct = compLambda(sz, a, true);
%   Wx_1 = dctn_fftw(L_dct.*dctn_fftw(x), 1);
%   Wx_2 = reshape(W_dct*x(:),sz);
%   disp( max( abs(Wx_1(:) - Wx_2(:)) ) )
%
%   W_fft = createQ(sz, 1, a, [], false);
%   L_fft = compLambda(sz, a, false);
%   Wx_1 = real(ifftn( L_fft.*fftn(x) ));
%   Wx_2 = reshape( W_fft*x(:), sz);
%   disp( max( abs(Wx_1(:) - Wx_2(:)) ) )

% $Id: compLambda.m 3724 2014-07-17 16:07:09Z johanl $

if nargin<2 || isempty(a), a=1; end
if nargin<3 || isempty(Neumann), Neumann=true; end

%number of dimensions
d = length(sz);

if ~isscalar(a) && length(a)~=d
  error('a must be scalar or have length(sz) elements')
end
if any(a<0), error('a must be non-negative'); end

%repeat a to match length of sz.
if isscalar(a), a=repmat(a, [1 length(sz)]); end

% should we use pi/sz (DCT) or 2*pi/sz (FFT)
if Neumann, pi_sz = pi./sz; else pi_sz = 2*pi./sz;  end

% Lambda tensor
Lambda = zeros(sz);
if d==1
  Lambda = a(1) * cos(pi_sz(1)*(0:(sz(1)-1))');
else
  for i = 1:d
    sz0 = ones(1,d);
    sz0(i) = sz(i);
    Lambda = bsxfun(@plus, Lambda, ...
                    a(i) * cos(pi_sz(i)*reshape(0:(sz(i)-1),sz0)));
  end
end
Lambda = 2*(sum(a)-Lambda);
