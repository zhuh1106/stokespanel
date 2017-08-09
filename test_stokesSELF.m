% test Stokes self evaluation
close all
clear all
v = 1;
side = 'i'; % test interior or exterior
lptype = 'd'; % test SLP or DLP
qntype = 'C'; % quadrature nodes, test gauss or chebyshev  
N = 100;

% set up source and target
% source: starfish domain
a = .3; w = 5;           % smooth wobbly radial shape params...
R = @(t) (1 + a*cos(w*t))*1; Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
% inside = @(z) abs(z)<R(angle(z));   % Boolean true if z inside domain
% inside = @(z) abs(z)<1.3;   % Boolean true if z inside domain
% outside = @(z) abs(z)>R(angle(z));   % Boolean true if z outside domain
s = quadr(s, N);

[A,Ag] = stokesselfevalm(s, N, lptype, side, qntype);
diff = Ag-A;
figure()
surf(diff(1:16:end,1:16:end));
