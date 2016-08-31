function Acorr = stokespanelcorm(px, rx, a, b, lptype, side, qntype)
% STOKESPANELCOR - Panel Correction Scheme for Stokes Equation
%
% Acorr = stokespanelcor(px, rx, a, b, lptype, side) gives special
% quadrature value matrix for panel correction.
% Inputs: px = target node, with or without struct both will work.
%         rx = source node, with or without struct both will work.
%         a = panel start, b = panel end, in complex plane.
%         lptype = SLP, 's', or DLP, 'd'.
%         side = interior, 'i', or exterior, 'e'.
% Output: Acorr is special quadrature value at target node px.
% Efficient only if multiple targs, since O(p^3).
% See Helsing-Ojala 2008 (special quadr Sec 5.1-2), Helsing 2009 mixed (p=16),
% and Helsing's tutorial demo11b.m LGIcompRecFAS()
% Hai 08/28/16

be = 2;     % factor by which to incr panel nodes for close eval
p = []; if isstruct(px), p.x = px.x; else p.x = px; end     % form target nodes struct
r = []; if isstruct(rx), r.x = rx.x; else r.x = rx; end     % form source nodes struct
% get struct for r with geometry info
if ~isfield(r,'nx'), r = quadr_panf(r, 1, qntype); end   
rf = quadr_panf(r, be, qntype); % struct with geometry info at fine nodes for close eval
num = numel(r.x);
Imn = interpmat(num, num*be, qntype);     % interpolation matrix
ta1M = [eye(num),zeros(num)];   ta2M = [zeros(num),eye(num)];   % matrix maps ta to ta1 and ta2

if lptype=='s'
    [A, A1, A2] = Sspecialquad(p,rf,a,b,side); % close-eval
    Acorr = A*Imn*( ta1M + 1i*ta2M)/2 ...
        + (A1+1i*A2)*Imn*( diag(real(r.x))*ta1M + diag(imag(r.x))*ta2M)/2 ...
        - diag(real(p.x))*( (A1+1i*A2)*Imn*ta1M)/2 ...
        - diag(imag(p.x))*( (A1+1i*A2)*Imn*ta2M)/2;
    % rewrite Stokes single layer potential as sum of Laplace potentials
else
    [A, A1, A2] = Dspecialquad(p,rf,a,b,side);  % close-eval
    ta11M = diag(real(r.nx)./r.nx)*( ta1M + 1i*ta2M);   % matrix maps ta to sigma (complex density function)
    ta12M = diag(imag(r.nx)./r.nx)*( ta1M + 1i*ta2M);
    Acorr = real(A*Imn*ta11M)+1i*real(A*Imn*ta12M) ...
        + (A1+1i*A2)*Imn*( diag(real(r.x))*ta1M + diag(imag(r.x))*ta2M)...
        - diag(real(p.x))*( (A1+1i*A2)*Imn*ta1M)...
        - diag(imag(p.x))*( (A1+1i*A2)*Imn*ta2M);
    % rewrite Stokes single layer potential as sum of Laplace potentials

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
Nn =  numel(find(inr));
if Nn ~= 0  % Criterion added by Hai Zhu 08/24/16 to ensure inr not empty
    for k=1:p
        P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); 
    end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf =  numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf ~= 0  % Criterion added by Hai Zhu 08/24/16 to ensure ifr not empty
    for k=2:p+1,
        P(k,ifr) = sum(repmat(wxp.*y.^(k-1),[1 Nf])./...
            (repmat(y,[1 Nf])-repmat(x(ifr).',[p 1])),1);
    end
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

function sf = quadr_panf(s, be, qntype)  
% set up quadrature on a closed segment
% QUADR_panf - set up quadrature (either coarse or fine nodes) on a segment struct
%
% sf = quadr_panf(s, be) gives quadrature on coarse or fine nodes.
% Inputs: s  = segment struct containing parametrization
%         be = factor by which to increase panel nodes
%  
% Outputs: sf - segment struct on fine node
% Hai 08/23/16

if ~isfield(s,'p')
    s.p=16; 
end
p = s.p; % default panel order
sf=[]; sf.p=be*s.p; pf=sf.p;
Imn = interpmat(p, pf, qntype);
sf.x = Imn*s.x;
if qntype=='G', [~, w, D] = gauss(pf); else [~, w, D] = cheby(pf); end 

sf.xp = D*sf.x; % velocities Z'(sf.x)
sf.xpp = D*sf.xp;   % acceleration Z''(sf.x)
sf.w = w;
sf.sp = abs(sf.xp); sf.tang = sf.xp./sf.sp; sf.nx = -1i*sf.tang;    % outward unit normals
sf.cur = -real(conj(sf.xpp).*sf.nx)./sf.sp.^2;
sf.ws = sf.w.*sf.sp; % speed weights
sf.wxp = sf.w.*sf.xp; % complex speed weights (Helsing's wzp)


function P = interpmat(n,m, qntype) % interpolation matrix from n-pt to m-pt Gauss nodes
% INTERPMAT - create interpolation matrix from n-pt to m-pt Gauss nodes
%
% P = interpmat(n,m) returns a m*n matrix which maps func values on n-pt Gauss-
% Legendre nodes on [-1,1] to values on m-pt nodes.
% Does it the Helsing way via backwards-stable ill-cond Vandermonde solve.
if m==n, P = eye(n); return, end
if qntype=='G', x = gauss(n); y = gauss(m); 
else x = cheby(n); y = cheby(m); end 
V = ones(n); for j=2:n, V(:,j) = V(:,j-1).*x; end % Vandermonde, original nodes
R = ones(m,n); for j=2:n, R(:,j) = R(:,j-1).*y; end % monomial eval matrix @ y
P = (V'\R')';                                       % backwards-stable solve

