function [ss,N] = triangle_test(x1, x2, x3, t, k, r)
% function to generate analytic expression of a smooth triangle
%
% Input:
% x1, x2, x3: corners of a triangle
% t:          each side of a triangle after smoothing is going to be a
%             piecewise defined function with 3 pieces, with starting and ending info
%             in parametrization stored in t.
% N:          number of quadrature points
% r:          radius
%
% Output:
% ss:         a structure contains all info of smoothed triangle
% N:          number of quadrature points

qtype = 'p';
qntype = 'C';
bdType = 'triangle';
idNum = 'middle';
N = k*16*9;

Zt1 = cInfBoundary2([x1; x2; x3], r, [t(1),t(2),t(3),t(4)], bdType, idNum);
Zt2 = cInfBoundary2([x2; x3; x1], r, [t(4),t(5),t(6),t(7)], bdType, idNum);
Zt3 = cInfBoundary2([x3; x1; x2], r, [t(7),t(8),t(9),t(10)], bdType, idNum);


%% assembly
ss.Z = @(t) (t<=2*pi/3) .* Zt1(t) + ((t>2*pi/3).* (t<=4*pi/3)) .* Zt2(t)...
        + (t>4*pi/3) .* Zt3(t);

%% quadrature pts get  
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
% figure()
% plot(real(ss.x),imag(ss.x))
% daspect([1 1 1])
% hold on


% x1 = [0,0];
% x2 = [1,1];
% % x3 = [1,1];
% x3 = [1,0];

% N = 128*4;
% t1 = 0; t2 = 2*pi/3; t3 = 4*pi/3; t4 = 2*pi;
