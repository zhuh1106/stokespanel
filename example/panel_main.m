function panel_main()

close all
clear all

xLim = 6;
r = 0.1;
nc = 2;
k = 2;
nps = k*[4;5;5;8;5;5;4;5;5;8;5;5];   % number of panels along each side
%% get geometry
sc = panel_confinedgeo_gen(xLim,nps,r,nc);
N = numel(sc.x);

% plot(real(sc.x),imag(sc.x),'.');
% hold on
% plot(real(sc.x(1)),imag(sc.x(1)),'*')
% axis equal
% 
%% set up rhs (boundary condition) and solve for density
f = [(1-imag(sc.x).^2).*(abs(imag(sc.x))<1);zeros(N,1)];
% quiver(real(sc.x),imag(sc.x), f(1:end/2),f(end/2+1:end), 0.5);


%% set up islands
nps = k*[8;3;8;3];   % number of panels along each side
si = {};
shift = -0.6; r = 0.1;
si{1} = panel_islands_gen(nps,r,shift,nc);
si{1}.n = 16*sum(nps);

%% set up tank shape
nps = k*[4;5;4;4;3;4;3;4;3;4;4;4];   % number of panels along each side
r = 0.08;
shift = [-0.3,0.4];
st = panel_tank_gen(nps,r,shift,nc);    % number of panels at corner   
st.n = 16*sum(nps);


%% set up target point
nx = 40; gx = ((1:nx)/nx*2-1)*xLim; ny = 40; gy = ((1:ny)/ny*2-1)*3;
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[INc, ONc] = inpolygon(real(zz),imag(zz),real(sc.x),imag(sc.x));
iic = INc & ~ONc;

i = 1;
[INi, ONi] = inpolygon(real(zz),imag(zz),real(si{i}.x),imag(si{i}.x));
iii = ~INi;
ii = iic & iii;

[INi, ONi] = inpolygon(real(zz),imag(zz),real(st.x),imag(st.x));
iii = ~INi;
ii = ii & iii;
t.x = zz(ii(:)); 

%% set up self evaluation system with confined geometry
[Adlp,~] = stokesselfevalm(sc,sc, N, 'd', 'i', 'C'); 
Ascsc = -eye(size(Adlp))/2 + Adlp;

%% set up self evaluation system with islands
[Aslp,~] = stokesselfevalm(si{i},si{i}, si{i}.n, 's', 'e', 'C');
[Adlp,~] = stokesselfevalm(si{i}, si{i}, si{i}.n, 'd', 'e', 'C');
Asisi = Aslp + Adlp + eye(size(Adlp))/2;

%% set up self evaluation system with tank
[Aslp,~] = stokesselfevalm(st,st, st.n, 's', 'e', 'C');
[Adlp,~] = stokesselfevalm(st, st, st.n, 'd', 'e', 'C');
Astst = Aslp + Adlp + eye(size(Adlp))/2;


%% from island to confined geometry
Aslp = stokescloseevalm(sc, si{i}, si{i}.n, 's', 'e', 'C');
Adlp = stokescloseevalm(sc, si{i}, si{i}.n, 'd', 'e', 'C');
Ascsi = Aslp + Adlp;

%% from tank to confined geometry
Aslp = stokescloseevalm(sc, st, st.n, 's', 'e', 'C');
Adlp = stokescloseevalm(sc, st, st.n, 'd', 'e', 'C');
Ascst = Aslp + Adlp;

%% from confined geometry to island
Asisc = stokescloseevalm(si{i}, sc, N, 'd', 'i', 'C');

%% from tank to island
Aslp = stokescloseevalm(si{i}, st, st.n, 's', 'e', 'C');
Adlp = stokescloseevalm(si{i}, st, st.n, 'd', 'e', 'C');
Asist = Aslp + Adlp;

%% from confined geometry to tank
Astsc = stokescloseevalm(st, sc, N, 'd', 'i', 'C');

%% from island to tank
Aslp = stokescloseevalm(st, si{i}, si{i}.n, 's', 'e', 'C');
Adlp = stokescloseevalm(st, si{i}, si{i}.n, 'd', 'e', 'C');
Astsi = Aslp + Adlp;


%% form system with rhs
A = [Ascsc,Ascsi,Ascst;Asisc,Asisi,Asist;Astsc,Astsi,Astst];
rhs = [f;zeros(2*si{i}.n,1);zeros(2*st.n,1)];
tau = A\rhs;
tauc = tau(1:2*N); taui = tau(2*N+1:2*(N+si{i}.n)); taut = tau(2*(N+si{i}.n)+1:end);


%
% figure()
% plot([sc.t;sc.t+2*pi;si{i}.t+4*pi;si{i}.t+6*pi],tau)
% 
u = nan*(1+1i)*zz;
u(ii) = stokescloseeval(t, sc, tauc, N, 'd', 'i', 'C');
u(ii) = u(ii) + stokescloseeval(t, si{i}, taui, si{i}.n, 's', 'e', 'C')...
              + stokescloseeval(t, si{i}, taui, si{i}.n, 'd', 'e', 'C')...
              + stokescloseeval(t, st, taut, st.n, 's', 'e', 'C')...
              + stokescloseeval(t, st, taut, st.n, 'd', 'e', 'C');


figure()
plot(real(sc.x),imag(sc.x),'.r'); hold on
plot(real(si{i}.x),imag(si{i}.x),'.r');
u0 = 1; u1c = min(max(real(u),-u0),u0); u2c = min(max(imag(u),-u0),u0);
quiver(xx,yy, u1c.*ii,u2c.*ii, 0.8);
% quiver(xx,yy, real(u).*ii,imag(u).*ii, 0.8);
axis equal

plot(real(st.x),imag(st.x),'.r');





