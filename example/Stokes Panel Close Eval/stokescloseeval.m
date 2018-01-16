function u = stokescloseeval(t, s, tau, N, lptype, side, qntype, normal_opt)

if nargin<5, lptype = 'd'; end     % 's' SLP or 'd' DLP test
if nargin<6, side = 'i'; end     % 'i'=interior or 'e'=exterior
if nargin<7, qntype = 'C'; end  % 'C'=chebyshev or 'G'=gauss
if nargin<8, normal_opt = 0; end
qtype = 'p'; s.p = 16;   % 'g'=global, 'p'=panels. s.p= panel order

%% set up panel-wise quadrature

if nargin<4, N = 300; end
% [s, N, np] = quadr(s,N,qtype,qntype); % set up bdry nodes (note outwards normal)
s_test = s;
[s, N, np] = quadr_pan(s, N, qtype, qntype) ;
np = numel(s.x)/s.p;
if normal_opt == 1
    s.nx = -s.nx; s.cur = -s.cur;
end
% sf = s; be = 2.0;        % factor by which to incr panel nodes for close eval
% sf.p = be*s.p; sf = quadr_pan(sf,be*N,qtype,qntype); % the fine nodes for close eval


%% Do native evaluation
if isempty(tau)
    disp('tau not defined!')
    % tau = [sin(2*pi - s.t);cos(2*pi - s.t)];
    tau = [sin(s.t);cos(s.t)]; % default tau for debug
end

if lptype=='s'
    Ag = SLPmatrix(t,s);
else
    Ag = DLPmatrix(t,s); 
end

% max(abs(test-test_mat*tau_test(1:64)))
% max(max(abs(Ag-test_mat)))
% max(abs(tau-tau_test(1:64)))

u_temp = Ag * tau; % native u eval grid (could be via FMM)
u_return = u_temp;
u = u_temp(1:end/2) + 1i*u_temp(end/2+1:end);


%% For each panel, replace close-native-eval by close-special-eval
if qtype == 'p'
tic, nst = 0; % eval on grid w/ close-eval scheme corrections

panlen = zeros(np); % lengths of panel (keep for later)
for k=1:np     % outer loop over panels
    j = (k-1)*s.p+(1:s.p);          % indices of nodes for kth panel
    panlen(k) = sum(s.ws(j));       % lengths of panel (keep for later)
    ik = (abs(t.x - (s.xlo(k)+s.xhi(k))/2) < 0.8*panlen(k) ...
        | abs(t.x - s.x(j(ceil(s.p/2)))) < panlen(k));
%         & abs(t.x - s.x(j(ceil(s.p/2))))<.5;
    
    % approx criterion for near-field of panel
    p.x = t.x(ik(:)); nst = nst + numel(p.x);  % col vec of targs near kth panel
    r.x = s.x(j); r.nx = s.nx(j); r.sp = s.sp(j); r.ws = s.ws(j); r.wxp = s.wxp(j);% r = this panel
    
    if numel(p.x)*numel(r.x)~=0
        % cancel native evaluation from panel k
        ta = tau([j,j+N]);                    % panel src nodes, density  

        if lptype=='s' 
            A =  SLPmatrix(p,r); 
        else
            A =  DLPmatrix(p,r);
        end
        u_temp = A*ta;
        u(ik) = u(ik) - (u_temp(1:end/2)+1i*u_temp(end/2+1:end));
        
%        jf = (k-1)*sf.p+(1:sf.p);       % indices of this fine panel
%        rf = []; rf.x = sf.x(jf); rf.nx = sf.nx(jf); rf.wxp = sf.wxp(jf); % fine panel
%        p.x = p.x(1);  % for debug dimension 1 case
        % add Helsing value (evalu`ated at fine nodes)
%        I = stokespanelcor( p.x, r.x, s.xlo(k), s.xhi(k), ta, lptype, side);   % density function as input, return correction value
        
        Acorr = stokespanelcorm( p.x, r.x, s.xlo(k), s.xhi(k), lptype, side, qntype);   % no density function, return correction matrix
        I1 = Acorr*ta;
        
        u(ik) = u(ik) + I1;     
        
    end

end
% uu = u_return(1:end/2) + 1i*u_return(end/2+1:end);
% abs(uu-u)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s, N, np] = quadr(s, N, qtype, qntype)  % set up quadrature on a closed segment
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

function A = DLPmatrix(t,s) % stokes double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;

if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
    t1 = real(s.tang); t2 = imag(s.tang);
    A(diagind(A)) = [t1.^2; t2.^2].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
    A(diagind(A)) = [t1.*t2; t2.*t1].*repmat(-s.cur,2,1)/2/pi;
    A = A(:,[1+end/2:end,1:end/2]);
end           % self? diagonal term for Stokes

A = A.*repmat(s.ws(:)', [2*M 2]);


function A = SLPmatrix(t,s) % double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
r = abs(d);    % dist matrix R^{MxN}

log_part = kron(eye(2),-log(r));

d1 = real(d)./r;
d2 = imag(d)./r;

cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

A = (log_part + cross_part).*repmat(s.ws(:)', [2*M 2])/4/pi;




function i = diagind(A) % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);


