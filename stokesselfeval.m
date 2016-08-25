function u = stokesselfeval(s, tau, N, lptype, side)

if nargin<4, lptype = 'd'; end     % 's' SLP or 'd' DLP test
if nargin<5, side = 'i'; end     % 'i'=interior or 'e'=exterior
qtype = 'p'; s.p = 16;   % 'g'=global, 'p'=panels. s.p= panel order

%% set up panel-wise quadrature

if nargin<3, N = 300; end
[s, N, np] = quadr(s,N,qtype); % set up bdry nodes (note outwards normal)

sf = s; be = 2.0;        % factor by which to incr panel nodes for close eval
sf.p = be*s.p; sf = quadr(sf,be*N,qtype); % the fine nodes for close eval


%% Do native evaluation
if nargin<2 || isempty(tau)
    disp('tau not defined!')
    % tau = [sin(2*pi - s.t);cos(2*pi - s.t)];
    tau = [sin(s.t);cos(s.t)]; % default tau for debug
end

if lptype=='s', Ag = SLPmatrix(s,s);
else Ag = DLPmatrix(s,s); end

u_temp = Ag * tau; % native u eval grid (could be via FMM)
u = u_temp(1:end/2) + 1i*u_temp(end/2+1:end);


%% For each panel, replace close-native-eval by close-special-eval
if qtype == 'p'
tic, nst = 0; % eval on grid w/ close-eval scheme corrections
Imn = interpmat(s.p, sf.p); % coarse-to-fine interp matrix, same for all panels
panlen = zeros(1,np); % lengths of panel (keep for later)
for k=1:np     % outer loop over panels
    j = (k-1)*s.p+(1:s.p);          % indices of nodes for kth panel
    panlen(k) = sum(s.ws(j));       % lengths of panel (keep for later)
    ik = (abs(s.x - (s.xlo(k)+s.xhi(k))/2) < 0.7*panlen(k) ...
        | abs(s.x - s.x(j(ceil(s.p/2)))) < 0.9*panlen(k));
    % approx criterion for near-field of panel
    p.x = s.x(ik(:)); nst = nst + numel(p.x);  % col vec of targs near kth panel
    r.x = s.x(j); r.nx = s.nx(j); r.sp = s.sp(j); r.ws = s.ws(j); % r = this panel
    if numel(p.x)*numel(r.x)~=0
        ta = tau([j,j+N]);                    % panel src nodes, density
        if lptype=='s', A =  SLPmatrix(p,r); else A =  DLPmatrix(p,r);end
        u_temp = A*ta;
        u(ik) = u(ik) - (u_temp(1:end/2)+1i*u_temp(end/2+1:end));      % cancel off local direct interaction

        I = stokespanelcor( p.x, r.x, s.xlo(k), s.xhi(k), ta, lptype, side);      
        u(ik) = u(ik) + I;        % add Helsing value
        
    end
end
end
u = [real(u);imag(u)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, N, np] = quadr(s, N, qtype)  % set up quadrature on a closed segment
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
    [~, ~, D] = gauss(N);
elseif qtype=='p'
    if ~isfield(s,'p'), s.p=16; end, p = s.p; % default panel order
    np = ceil(N/p); N = p*np;      % np = # panels
    s.tlo = (0:np-1)'/np*2*pi; s.xlo = s.Z(s.tlo); % panel start params, locs
    s.thi = (1:np)'/np*2*pi; s.xhi = s.Z(s.thi);  % panel end params, locs
    pt = 2*pi/np;                  % panel size in parameter
    t = zeros(N,1); s.w = t;
    [x, w, D] = gauss(p);
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

function A = DLPmatrix(t,s) % stokes double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
r2 = abs(d).^2;    % dist^2 matrix R^{MxN}

dot_part = real(ny./d)./r2/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dot_part(isinf(dot_part)) = 0;%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;

% if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
%     t1 = real(s.tang); t2 = imag(s.tang);
%     A(diagind(A)) = [t1.^2; t2.^2].*repmat(-s.cur,2,1)/2/pi;
%     A = A(:,[1+end/2:end,1:end/2]);
%     A(diagind(A)) = [t1.*t2; t2.*t1].*repmat(-s.cur,2,1)/2/pi;
%     A = A(:,[1+end/2:end,1:end/2]);
% end           % self? diagonal term for Stokes

A = A.*repmat(s.ws(:)', [2*M 2]);


function A = SLPmatrix(t,s) % double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
r = abs(d);    % dist matrix R^{MxN}
%%%%%%%%%%%%%
r(r==0) = 1;%
%%%%%%%%%%%%%
log_part = kron(eye(2),-log(r));

d1 = real(d)./r;
d2 = imag(d)./r;

cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

A = (log_part + cross_part).*repmat(s.ws(:)', [2*M 2])/4/pi;

function [A, A1, A2] = Sspecialquad(t,s,a,b,side)
% SSPECIALQUAD - SLP val+grad close-eval Helsing "special quadrature" matrix
%
% [A] = Sspecialquad(t,s,a,b) returns
% [A An] = Sspecialquad(t,s,a,b) also gives target normal-derivs (needs t.nx)
% [A A1 A2] = Sspecialquad(t,s,a,b) also gives target x,y-derivs
% Inputs: t = target seg struct (with column-vec t.x targets in complex plane)
%         s = src node seg struct (with s.x, s.w, s.nx)
%         a = panel start, b = panel end, in complex plane.
% Output: A (n_targ * n_src) is source-to-target value matrix
%         An or A1, A2 = source to target normal-deriv (or x,y-deriv) matrices
% Efficient only if multiple targs, since O(p^3).
% See Helsing-Ojala 2008 (special quadr Sec 5.1-2), Helsing 2009 mixed (p=16),
% and Helsing's tutorial demo11b.m LGIcompRecFAS()
if nargin<5, side = 'i'; end     % interior or exterior
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p+1,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
% compute P up to p+1 instead of p as in DLP, since q_k needs them:
if side == 'i'
    P(1,:) = 1i*pi/2 + log((1-x)./(1i*(-1-x)));  % initialize p_1 for all targs int
elseif side == 'e'
    P(1,:) = -1i*pi/2 + log((1-x)./(-1i*(-1-x)));  % initialize p_1 for all targs ext
end
% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
for k=1:p, P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); end  % recursion for p_k
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf =  numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
for k=2:p+1,
    P(k,ifr) = sum(repmat(wxp.*y.^(k-1),[1 Nf])./...
        (repmat(y,[1 Nf])-repmat(x(ifr).',[p 1])),1);
end
% note this is slower than recurrence
Q = zeros(p,N); % compute q_k from p_k via Helsing 2009 eqn (18)... (p even!)
if side == 'i'
    Q(1:2:end,:) = P(2:2:end,:) - repmat(log(1i*(1-x.'))+log(1i*(-1-x.')),[p/2 1]);
    % (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
    Q(2:2:end,:) = P(3:2:end,:) - repmat(log(1i*(1-x.'))-log(1i*(-1-x.')),[p/2 1]);
elseif side == 'e'
    Q(1:2:end,:) = P(2:2:end,:) - repmat(log(-1i*(1-x.'))+log(-1i*(-1-x.')),[p/2 1]);
    % (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
    Q(2:2:end,:) = P(3:2:end,:) - repmat(log(-1i*(1-x.'))-log(-1i*(-1-x.')),[p/2 1]);
end
Q = Q.*repmat(1./(1:p)',[1 N]); % k even
warning('off','MATLAB:nearlySingularMatrix'); % solve for special weights...
A = real((V.'\Q).'.*repmat((1i*s.nx)',[N 1])*zsc)/(2*pi*abs(zsc));
A = A*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N 1]); % unscale, yuk
if nargout>1
    P = P(1:end-1,:);  % trim P back to p rows since kernel is like DLP
    Az = (V.'\P).'*(1/(2*pi)).*repmat((1i*s.nx)',[N 1]); % solve spec wei
    A1 = real(Az); A2 = -imag(Az);
end

function [A, A1, A2] = Dspecialquad(t,s,a,b,side)
% DSPECIALQUAD - DLP val+grad close-eval Helsing "special quadrature" matrix
%
% [A] = Dspecialquad(t,s,a,b)
% [A An] = Dspecialquad(t,s,a,b) also gives target normal-derivs (needs t.nx)
% [A A1 A2] = Dspecialquad(t,s,a,b) also gives target x,y-derivs
% Inputs: t = target seg struct (with column-vec t.x targets in complex plane)
%         s = src node seg struct (with s.x, s.w; amazingly, s.nx not used!)
%         a = panel start, b = panel end, in complex plane.
% Output: A (n_targ * n_src) is source-to-target value matrix
%         An or A1, A2 = source to target normal-deriv (or x,y-deriv) matrices
% Efficient only if multiple targs, since O(p^3).
% See Helsing-Ojala 2008 (special quadr Sec 5.1-2), Helsing 2009 mixed (p=16),
% and Helsing's tutorial demo11b.m M1IcompRecFS()
if nargin<5, side = 'i'; end     % interior or exterior
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
%figure; plot(x,'.'); hold on; plot(y,'+-'); plot([-1 1],[0 0],'ro'); % debug
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if N*p==0
    A = 0; A1=0; A2=0;
    return
end
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
if side == 'i'
    P(1,:) = 1i*pi/2 + log((1-x)./(1i*(-1-x)));  % initialize p_1 for all targs int
elseif side == 'e'
    P(1,:) = -1i*pi/2 + log((1-x)./(-1i*(-1-x)));  % initialize p_1 for all targs ext
end
% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
if N>1 || (N==1 && inr==1) % Criterion added by Bowei Wu 03/05/15 to ensure inr not empty
    for k=1:p-1, P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf = numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf>0 % Criterion added by Bowei Wu 03/05/15 to ensure ifr not empty
    for k=2:p, P(k,ifr) = sum(repmat(wxp.*y.^(k-1),[1 Nf])./...
            (repmat(y,[1 Nf])-repmat(x(ifr).',[p 1])),1); end   % note this is slower than recurrence
end
warning('off','MATLAB:nearlySingularMatrix');
% A = real((V.'\P).'*(1i/(2*pi)));         % solve for special quadr weights
A = ((V.'\P).'*(1i/(2*pi)));         % do not take real for the eval of Stokes DLP non-laplace term. Bowei 10/19/14
if nargout>1
    R =  -(kron(ones(p,1),1./(1-x.')) + kron((-1).^(0:p-1).',1./(1+x.'))) +...
        repmat((0:p-1)',[1 N]).*[zeros(1,N); P(1:p-1,:)];  % hypersingular kernel weights of Helsing 2009 eqn (14)
    Az = (V.'\R).'*(1i/(2*pi*zsc));  % solve for targ complex-deriv mat & rescale
    %if nargout==2
    %An = real(repmat(t.nx, [1 p]) .* Az);  % to check - rotate nx ***
    %else
    A1 = real(Az); A2 = -imag(Az);   % note sign for y-deriv from C-deriv
    %end
end

function i = diagind(A)
% function i = diagind(A)
%
% return diagonal indices of a square matrix, useful for changing a diagonal
% in O(N) effort, rather than O(N^2) if add a matrix to A using matlab diag()
%
% barnett 2/6/08
N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);

function P = interpmat(n,m) % interpolation matrix from n-pt to m-pt Gauss nodes
% INTERPMAT - create interpolation matrix from n-pt to m-pt Gauss nodes
%
% P = interpmat(n,m) returns a m*n matrix which maps func values on n-pt Gauss-
% Legendre nodes on [-1,1] to values on m-pt nodes.
% Does it the Helsing way via backwards-stable ill-cond Vandermonde solve.
if m==n, P = eye(n); return, end
x = gauss(n); y = gauss(m);
V = ones(n); for j=2:n, V(:,j) = V(:,j-1).*x; end % Vandermonde, original nodes
R = ones(m,n); for j=2:n, R(:,j) = R(:,j-1).*y; end % monomial eval matrix @ y
P = (V'\R')';                                       % backwards-stable solve

