function test_stokesBVP
% TEST_STOKESBVP - Test Stokes Close Evaluation by Solving a BVP
%
%        test_stokesBVP returns log scale err plot for solution obtained by
%        using close evaluation against exact solution
% Variables : t = target node( subset of nodes of s)
%             s = source node( boundary), with or without struct both will work.
%             lptype = SLP, 's', or DLP, 'd'.
%             side = interior, 'i', or exterior, 'e'.
%             qntype = chebyshev, 'C', or gauss, 'G'
% Output:   err plot
% Hai 08/26/17

% test Stokes BVP
v = 1;
side = 'i'; % test interior or exterior
lptype = 'd'; % test SLP or DLP
qntype = 'C'; % quadrature nodes, test gauss or chebyshev  
N = 600;

% set up source and target
% source: starfish domain
a = .3; w = 5;           % smooth wobbly radial shape params...
R = @(t) (1 + a*cos(w*t))*1; Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s = quadr(s, N);

% target
nx = 150; gx = ((1:nx)/nx*2-1)*1.5; ny = 150; gy = ((1:ny)/ny*2-1)*1.5; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);


% generate the exact solution

% set up locations of point forces
alpha = linspace(-pi/4,pi/4,5);
y_force = [];                             
if side == 'e'
    y_force.x = .3*exp(1i*alpha).';   % location of the point forces (exterior)
else
    y_force.x = 2*exp(1i*alpha).'; % location of the point forces (interior)
end
% pt_force = [[1;1;0;1;0];[0;0;1;0;1]];  % magnitude & direction of the point forces
pt_force = [[1;1;0;1;0];[1;0;-1;0;1]];  % magnitude & direction of the point forces (sample 2)


% use stokeslet to eval exact soln
fhom = nan*(1+1i)*zz; % exact soln
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x));
if side == 'i'
    ii = IN & ~ON;
elseif side == 'e'
    ii = ~IN;
end
t.x = zz(ii(:));  % eval pts only on one side
A = stokesletmatrix(t,y_force);
f_temp = A*pt_force;
fhom(ii(:)) = f_temp(1:end/2) + 1i*f_temp(end/2+1:end); % the exact soln (complex form)

% plot the exact soluntion and the point forces
if v==1
    figure(3); streamslice(gx,gy,real(fhom),imag(fhom));
    hold on; plot(real(s.Z(s.t)),imag(s.Z(s.t)),'r');
    title('Exact Soln')
    plot(y_force.x,'.','MarkerSize',10,'LineWidth',10)
    quiver(real(y_force.x),imag(y_force.x),.2*pt_force(1:end/2),.2*pt_force(end/2+1:end))
    axis equal tight
end

% 2 convergence test

Nn = 10;
for NN = 1:Nn
    N = 100*NN;
    % solve the BVP
    [s, N, np] = quadr_pan(s,N,'p',qntype); % set up bdry nodes (note outwards normal)
    N  
    % get bdry condition
    f = stokesletmatrix(s,y_force)*pt_force; % Dirichlet bdry condition
    
    %% solve for tau using self close evaluation matrix
    warning('off','MATLAB:nearlySingularMatrix')
    if lptype == 's'
        if side == 'e'
            [As,~] = stokesselfevalm(s,s, N, 's', side, qntype);
            tau = As\f;
        elseif side == 'i'
            [As,~] = stokesselfevalm(s,s, N, 's', side, qntype);
            tau = As\f;
        end
    elseif lptype == 'd'
        if side == 'e'
            [Ad,~] = stokesselfevalm(s,s, N, lptype, side, qntype);
            [As,~] = stokesselfevalm(s,s, N, 's', 'e', qntype);
            tau = (eye(size(Ad))/2 + Ad + As)\f;
        elseif side == 'i'
            [A,~] = stokesselfevalm(s,s, N, lptype, side, qntype);
            tau = (-eye(size(A))/2 + A)\f;
        end
    end
    
    %% evaluate velocity inside computation domain using close evaluation scheme
    u = nan*(1+1i)*zz;
    if lptype == 's'
        u(ii) = stokescloseeval(t, s, tau, N,lptype,side,qntype);
        if side == 'i'
            % unique up to a rigid body motion
            z0 = [-.1+0i;(.2+0.1i)*exp(2i*pi/5)];  % two test pts inside
            u0 = stokescloseeval(struct('x',z0), s, tau, N,lptype,side,qntype); % vel @ test pts
            ue0 = stokesletmatrix(struct('x',z0),y_force)*pt_force; % exact soln
            ue0 = ue0(1:end/2) + 1i*ue0(end/2+1:end); % make complx rep @ test pts
            du0 = ue0-u0; uconst = du0(1); urot = imag((du0(2)-du0(1))/(z0(2)-z0(1))); % get 3 params
            u(ii) = u(ii) + uconst + 1i*urot*(t.x-z0(1)); % correct the evaluated vel field
        end
    elseif lptype == 'd'
        if side == 'e'
            u(ii) = stokescloseeval(t, s, tau, N,'s',side, qntype)+stokescloseeval(t, s, tau, N,'d',side, qntype);
        elseif side == 'i'
            u(ii) = stokescloseeval(t, s, tau, N, lptype, side, qntype);
        end
    end
    figure(1),clf,imagesc(gx,gy,log10(abs(u-fhom)))
    colorbar
    err(NN) = max(abs(u(:)-fhom(:)));
end

figure(1),clf, imagesc(gx,gy,log10(abs(u-fhom))), colorbar, title('log10 err in |u|'), axis equal tight
figure(2),clf, semilogy(100*(1:Nn),err,'o'), title('BVP conv')

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = stokesletmatrix(t,s) % stokeslet matrix
% t = target seg (x,nx cols), s = src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
r = abs(d);    % dist matrix R^{MxN}

log_part = kron(eye(2),-log(r));

d1 = real(d)./r;
d2 = imag(d)./r;

cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

A = (log_part + cross_part)/4/pi;

function A = stokesletpmatrix(t,s) % normal derive of stokeslet matrix
% t = target seg (x,nx cols), s = src seg
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(t.nx, [1 N]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = -real(ny./d)./r2/pi;
dot_part = repmat(dot_part,2,2);

d1 = real(d);
d2 = imag(d);
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];
A = dot_part.*cross_part;


function A = DLPmatrix(t,s) % double-layer kernel matrix
% t = target seg (x,nx cols), s = src seg
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

log_part = kron(eye(2),-log(r).*repmat(s.w(:)', [M 1]));

d1 = real(d)./r;
d2 = imag(d)./r;

cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2].*repmat(s.ws(:)', [2*M 2]);

A = (log_part + cross_part)/4/pi;

function A = SLPmatrixp(t,s) % normal deriv of single-layer kernel self matrix
% t = target seg (x,nx cols), s = src seg
% No jump included on self-interaction.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(t.nx, [1 N]);      % identical rows given by src normals
r2 = d.*conj(d);    % dist^2 matrix R^{MxN}

dot_part = -real(ny./d)./r2/pi;
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


function S = SLPselfmatrix(s) % Stokes SLP via Kress-split Nyst matrix
% Modified based on Laplace case (Barnett 10/18/13). Wu 10/13/14.

% log part
N = numel(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);  % C-# displacements mat, t-s
logpart = -log(abs(d)) + circulant(0.5*log(4*sin([0;s.t(1:end-1)]/2).^2)); % peri log
logpart(diagind(logpart)) = -log(s.sp);                       % diagonal limit
m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2; % Kress Rj(N/2)/4pi
logpart = logpart/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
logpart = logpart .* repmat(s.sp.',[N 1]);  % include speed factors (not 2pi/N weights)
logpart = logpart/2; % Stokes prefactor

% cross product part
r = sqrt(d.*conj(d));    % dist matrix R^{MxN}
d1 = real(d)./r;
d2 = imag(d)./r;
cross_part = [d1.^2, d1.*d2; d1.*d2, d2.^2];

t1 = real(s.tang); t2 = imag(s.tang);
cross_part(diagind(cross_part)) = [t1.^2; t2.^2];
cross_part = cross_part(:,[1+end/2:end,1:end/2]);
cross_part(diagind(cross_part)) = [t1.*t2; t2.*t1];
cross_part = cross_part(:,[1+end/2:end,1:end/2]);
    
cross_part = cross_part.*repmat(s.ws(:)', [2*N 2])/4/pi;

S = kron(eye(2),logpart)+cross_part;

function A = circulant(x)
% function A = circulant(x)
%
% return square circulant matrix with first row x
% barnett 2/5/08
x = x(:);
A = toeplitz([x(1); x(end:-1:2)], x);
