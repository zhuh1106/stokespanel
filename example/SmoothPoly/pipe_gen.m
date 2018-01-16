close all
clear all
r = 0.2;

xu = [0,1;1,0;2,0;3,1];
[su, N] = pipe_test(xu, r);

xl = [3,-2;2,-1;1,-1;0,-2];
[sl, N] = pipe_test(xl, r);

%% figure of cInf approx to a pipe
figure()
% plot(real(ss.x),imag(ss.x),'.')
plot(real(su.x),imag(su.x),'.')
hold on
plot(real(sl.x),imag(sl.x),'.')
daspect([1 1 1])

figure()
plot(linspace(0,2*pi,N),su.cur)
hold on
plot(linspace(0,2*pi,N),sl.cur)

% x = [0,3;1,1;3,2;2,-1;3,-3;0,-2;-3,-3;-2,-1;-3,2;-1,1];
% [ss,N] = polygon_test(x,r);
% ss;

%% upper and lower wall generated
% xu = [0,1;2,1;4,3;8,3;10,1;12,1];
xu = [12,1;10,1;8,3;4,3;2,1;0,1];
[su, N] = pipe_test(xu, r);

% xl = [12,-1;10,-1;8,-3;4,-3;2,-1;0,-1];
xl = [0,-1;2,-1;4,-3;8,-3;10,-1;12,-1];
[sl, N] = pipe_test(xl, r);


%% island inside channel generated

x = [2,0;9,1;8,2;4,2];
[ssu,N2] = polygon_test(x,r);

x = [3,-1;4,-2;8,-2;10,0];
[ssl,N2] = polygon_test(x,r);


%% plot for visualization
figure()
% plot(real(ss.x),imag(ss.x),'.')
plot(real(su.x),imag(su.x),'.')
hold on
plot(real(sl.x),imag(sl.x),'.')
plot(real(ssu.x),imag(ssu.x),'o')
plot(real(ssl.x),imag(ssl.x),'o')
daspect([1 1 1])

figure()
plot(linspace(0,2*pi,N),su.cur)
hold on
plot(linspace(0,2*pi,N),sl.cur)


%% set up source and target
% source: the tube bdry (orientation of source is such that targets are on the left)
su; sl; ssu; ssl;

% target
nx = 30; gx = (1:nx)/nx*12; ny = 15; gy = (1:ny)/ny*8-4; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);

[INp, ONp] = inpolygon(real(zz),imag(zz),[real(su.x);real(sl.x)],[imag(su.x);imag(sl.x)]);
INp;

[INi1, ONi1] = inpolygon(real(zz),imag(zz),real(ssu.x),imag(ssu.x));

[INi2, ONi2] = inpolygon(real(zz),imag(zz),real(ssl.x),imag(ssl.x));

[m,n] = size(INp);

% ii = reshape(INp,m*n,1) && reshape(~INi1,m*n,1) && reshape(~INi2,m*n,1);
% ii = reshape(ii,m,n);
ii = INp & (~INi1) & (~INi2);
ii;
t=[]; t.x = zz(ii(:));
% 
%% self convergence test
% N = 1500;
% f = nan*zz;
% lptype = 's'; % test SLP
% % lptype = 'd'; % test DLP
% tau = []; % tau need to be defined on the panel quadr nodes, here used default test
% f(ii) = stokescloseeval(t, su, tau, N,lptype);
% f(ii) = f(ii) + stokescloseeval(t, sl, tau, N,lptype);
% % 
% clear err
% k = 1;
% for N = 100:50:600
%     u = nan*zz;
%     u(ii) = stokescloseeval(t, su, tau, N,lptype);
%     u(ii) = u(ii) + stokescloseeval(t, sl, tau, N,lptype);
% %     imagesc(gx,gy,abs(u-f))
% %     colorbar
%     err(k) = max(abs(u(:)-f(:)));
%     k = k+1;
% end
% %% plot fig
% figure(1),imagesc(gx,gy,log10(abs(u-f))), colorbar, title('log10 err in |u|')
% figure(2), semilogy(100:50:600,err,'o'), title('self conv, against N=1500')
% 




