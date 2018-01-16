function ss = panel_confinedgeo_gen(xLim,nps,r,nc)
addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

%% set up source

% generate microfludic chip geometry
% xLim = 6; 
x1 = [-xLim,1]; x2 = [-xLim,-1]; x3 = [-4,-1]; x4 = [-3,-3];
x5 = [3,-3]; x6 = [4,-1]; x7 = [xLim,-1]; x8 = [xLim,1];
x9 = [4,1]; x10 = [3,3]; x11 = [-3,3]; x12 = [-4,1];
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];

% nps = [3;4;5;6;5;4;3;4;5;6;5;4];   % number of panels along each side
np = sum(nps); theta = 2*pi/np;
[polySize, ~] = size(x);
t = zeros(polySize,4);
t(1,:) = [0,nc*theta,nps(1)/np*2*pi-nc*theta,nps(1)/np*2*pi];
for i = 2:polySize
    t(i,:) = [t(i-1,end),t(i-1,end)+nc*theta,t(i-1,end)+(nps(i)-nc)*theta,t(i-1,end)+nps(i)*theta];
end

N = 16*np;
[ss,~] = polygon_test2(x,r,N,t);
%{
plot(real(ss.x),imag(ss.x),'.');
hold on
plot(real(ss.x(1)),imag(ss.x(1)),'*')
axis equal

%% set up rhs (boundary condition) and solve for density
f = [(1-imag(ss.x).^2).*(abs(imag(ss.x))<1);zeros(N,1)];
quiver(real(ss.x),imag(ss.x), f(1:end/2),f(end/2+1:end), 0.5);

%% set up target point
nx = 40; gx = ((1:nx)/nx*2-1)*xLim; ny = 40; gy = ((1:ny)/ny*2-1)*3;
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(ss.x),imag(ss.x));
ii = IN & ~ON;
t.x = zz(ii(:)); 

[A,~] = stokesselfevalm(ss,ss, N, 'd', 'i', 'C');
tau = (-eye(size(A))/2 + A)\f;
figure()
plot([ss.t;ss.t+2*pi],tau)

u = nan*(1+1i)*zz;
u(ii) = stokescloseeval(t, ss, tau, N, 'd', 'i', 'C');
figure()
plot(real(ss.x),imag(ss.x),'.r'); hold on
% u0 = 1; u1c = min(max(real(u),-u0),u0); u2c = min(max(imag(u),-u0),u0);
% quiver(xx,yy, u1c.*ii,u2c.*ii, 0.1);
quiver(xx,yy, real(u).*ii,imag(u).*ii, 0.5);
axis equal
%}
