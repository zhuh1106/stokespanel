%% set up source and target
% source: the tube bdry (orientation of source is such that targets are on the left)
epsilon = 0.1;
a=4;
b=1.5;
w = pi/2/6;
U = @(t) 2*pi-t + 1i*(3/8*w*tanh(epsilon*((pi-t).^a-(pi/b)^a))+7/8*w);
D = @(t) t - 1i*(3/8*w*tanh(epsilon*((t-pi).^a-(pi/b)^a))+7/8*w);
su.Z = U; % upper wall of the tube
sd.Z = D; % lower wall of the tube

% target
nx = 300; gx = (1:nx)/nx*2*pi; ny = 150; gy = (1:ny)/ny*2-1; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
inside = @(z) imag(z)<imag(U(real(z))) & imag(z)>imag(D(real(z)));   % Boolean true if z inside domain
ii = inside(zz);
t=[]; t.x = zz(ii(:));

%% self convergence test
N = 1500;
f = nan*zz;
lptype = 's'; % test SLP
% lptype = 'd'; % test DLP
tau = []; % tau need to be defined on the panel quadr nodes, here used default test
f(ii) = stokescloseeval(t, su, tau, N,lptype);
f(ii) = f(ii) + stokescloseeval(t, sd, tau, N,lptype);

clear err
k = 1;
for N = 100:50:600
    u = nan*zz;
    u(ii) = stokescloseeval(t, su, tau, N,lptype);
    u(ii) = u(ii) + stokescloseeval(t, sd, tau, N,lptype);
    err(k) = max(abs(u(:)-f(:)));
    k = k+1;
end
%% plot fig
figure(1),imagesc(gx,gy,log10(abs(u-f))), colorbar, title('log10 err in |u|')
figure(2), semilogy(100:50:600,err,'o'), title('self conv, against N=1500')


