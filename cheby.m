function [x, w, D] = cheby(N)
% based on chebyshev extreme cheb.m from spectral methods in matlab

% chebyshev nodes
theta = pi*(2*(1:N)'-1)/(2*N);
x = -cos(theta);
% chebyshev weights
l = floor(N/2)+1;
K = 0:N-l;   
v = [2*exp(1i*pi*K/N)./(1-4*K.^2)  zeros(1,l)];
w = real(ifft(v(1:N) + conj(v(N+1:-1:2))))';
% spectral differentiation matrix
X = repmat(x,1,N);
dX = X-X';
a = prod(dX+eye(N),2);
D = (a*(1./a)')./(dX+eye(N));
D = D - diag(sum(D,2));