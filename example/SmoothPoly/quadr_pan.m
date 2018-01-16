function [s, N, np] = quadr_pan(s, N, qtype, qntype)  % set up quadrature on a closed segment
% QUADR - set up quadrature (either global or panel-based) on a segment struct
%
% [s N] = quadr(s, N, qtype) adds quadrature input to segment struct s.
% Inputs:
%  s - segment struct containing parametrization
%  N - requested number of nodes
%  qtype - quadrature type: 'g' global periodic trapezoid rule
%                           'p' panel-wise Gauss-Legendre w/ s.p nodes per pan
% Outputs: s - amended segment struct
%          N - actual number of nodes (is within s.p of requested N)
%          np - number of panels
% Notes: 1) Note sign change in normal sense vs periodicdirpipe.m
% 2) Non-adaptive for now.  Barnett 6/10/13
if qtype=='g'
    t = (1:N)'/N*2*pi;
    s.tlo = 0; s.thi = 2*pi; s.p = N; s.w = 2*pi/N*ones(N,1); np=1; % 1 big panel
    s.xlo = s.Z(s.tlo); s.xhi = s.Z(s.thi);
    if qntype=='G', [~, ~, D] = gauss(N); else [~, ~, D] = cheby(N); end
    
elseif qtype=='p'
    if ~isfield(s,'p'), s.p=16; end, p = s.p; % default panel order
    np = ceil(N/p); N = p*np;      % np = # panels
    s.tlo = (0:np-1)'/np*2*pi; s.xlo = s.Z(s.tlo); % panel start params, locs
    s.thi = (1:np)'/np*2*pi; s.xhi = s.Z(s.thi);  % panel end params, locs
    pt = 2*pi/np;                  % panel size in parameter
    t = zeros(N,1); s.w = t;
    if qntype=='G', [x, w, D] = gauss(p); else [x, w, D] = cheby(p); end   
    D = D*2/pt;
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        t(ii) = s.tlo(i) + (1+x)/2*pt; s.w(ii) = w*pt/2; % nodes weights this panel
    end
end
s.x = s.Z(t);
s.xp = zeros(length(s.x),1);
s.xpp = zeros(length(s.x),1);

if isfield(s,'Zp'), s.xp = s.Zp(t);
else
    if qtype == 'p'
        for i=1:np
            ii = (i-1)*p+(1:p); % indices of this panel
            s.xp(ii) = D*s.x(ii);
        end
    else
        s.xp = D*s.x;
    end
end
if isfield(s,'Zpp'), s.xpp = s.Zpp(t);
else
    if qtype == 'p'
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        s.xpp(ii) = D*s.xp(ii);
    end
    else
        s.xpp = D*s.xp;
    end
end

s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang;
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2;
s.ws = s.w.*s.sp; % speed weights
s.t = t; s.wxp = s.w.*s.xp; % complex speed weights (Helsing's wzp)
s.w = s.ws;
s.cw = 1i*s.nx.*s.w;