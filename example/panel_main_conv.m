function [u,sc,si] = panel_main_conv(nc, k, npsc, npsi)

% close all

xLim = 6;
r = 0.1;

%% get geometry
nps = npsc*k;
sc = panel_confinedgeo_gen(xLim,nps,r,nc);
sc.n = numel(sc.x);
sc.np = sum(nps);
N = numel(sc.x);

%% set up rhs (boundary condition) and solve for density
f = [(1-imag(sc.x).^2).*(abs(imag(sc.x))<1);zeros(N,1)];

%% set up islands
option = 1;
si = islands(npsi,k,nc,option);

numOfI = 4;


%% set up target point
nx = 80; gx = ((1:nx)/nx*2-1)*xLim; ny = 40; gy = ((1:ny)/ny*2-1)*3;
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[INc, ONc] = inpolygon(real(zz),imag(zz),real(sc.x),imag(sc.x));
ii = INc & ~ONc;
numOfQ = sc.n;
for i = 1:numOfI
    [INi, ONi] = inpolygon(real(zz),imag(zz),real(si{i}.x),imag(si{i}.x));
    iii = ~INi;
    ii = ii & iii;
    numOfQ = numOfQ + si{i}.n;
end
t.x = zz(ii(:));


%% set up matrix
A = zeros(2*numOfQ);
incr_i = 0;
for i = 0:numOfI
    incr_j = 0;
    for j = 0:numOfI
        if (i==j)
            if (j==0)
                [Adlp,~] = stokesselfevalm(sc,sc, N, 'd', 'i', 'C'); 
                Asisj = -eye(size(Adlp))/2 + Adlp;
                A(incr_i+1:incr_i+2*sc.n, incr_j+1:incr_j+2*sc.n) = Asisj;
                incr_j = incr_j + 2*sc.n;
            else
                [Aslp,~] = stokesselfevalm(si{i}, si{i}, si{i}.n, 's', 'e', 'C');
                [Adlp,~] = stokesselfevalm(si{i}, si{i}, si{i}.n, 'd', 'e', 'C');
                Asisj = Aslp + Adlp + eye(size(Adlp))/2;
                A(incr_i+1:incr_i+2*si{i}.n, incr_j+1:incr_j+2*si{j}.n) = Asisj;
                incr_j = incr_j + 2*si{j}.n;
            end
        else
            if (j==0)
                Asisj = stokescloseevalm(si{i}, sc, N, 'd', 'i', 'C');
                A(incr_i+1:incr_i+2*si{i}.n, incr_j+1:incr_j+2*sc.n) = Asisj;
                incr_j = incr_j + 2*sc.n;
            elseif (i==0)
                Aslp = stokescloseevalm(sc, si{j}, si{j}.n, 's', 'e', 'C');
                Adlp = stokescloseevalm(sc, si{j}, si{j}.n, 'd', 'e', 'C');
                Asisj = Aslp + Adlp;
                A(incr_i+1:incr_i+2*sc.n, incr_j+1:incr_j+2*si{j}.n) = Asisj;
                incr_j = incr_j + 2*si{j}.n;
            else
                Aslp = stokescloseevalm(si{i}, si{j}, si{j}.n, 's', 'e', 'C');
                Adlp = stokescloseevalm(si{i}, si{j}, si{j}.n, 'd', 'e', 'C');
                Asisj = Aslp + Adlp;
                A(incr_i+1:incr_i+2*si{i}.n, incr_j+1:incr_j+2*si{j}.n) = Asisj;
                incr_j = incr_j + 2*si{j}.n;
            end
        end
    end
    if (i==0)
        incr_i = incr_i + 2*sc.n;
    else
        incr_i = incr_i + 2*si{i}.n;
    end
end

for i = 1:numOfI
    f = [f;zeros(2*si{i}.n,1)]; 
end
tau = A\f;
tauc = tau(1:2*N); 
incr = 2*N;
for i = 1:numOfI
    si{i}.tau = tau(incr+1:incr+2*si{i}.n);
    incr = incr + 2*si{i}.n;
end


u = nan*(1+1i)*zz;
u(ii) = stokescloseeval(t, sc, tauc, N, 'd', 'i', 'C');
for i = 1:numOfI
    u(ii) = u(ii) + stokescloseeval(t, si{i}, si{i}.tau, si{i}.n, 's', 'e', 'C')...
              + stokescloseeval(t, si{i}, si{i}.tau, si{i}.n, 'd', 'e', 'C');
end

figure()
plot(real(sc.x),imag(sc.x),'.r'); hold on
u0 = 1; u1c = min(max(real(u),-u0),u0); u2c = min(max(imag(u),-u0),u0);
quiver(xx,yy, u1c.*ii,u2c.*ii, 0.8);
% quiver(xx,yy, real(u).*ii,imag(u).*ii, 0.8);
axis equal
for i = 1:numOfI
    plot(real(si{i}.x),imag(si{i}.x),'.r');
end

function si = islands(npsi,k,nc,option)

switch option
    case 1
        nps = k*npsi;   % number of panels along each side
        si = {};
        shift = [-0.75,-0.75]; r = 0.1;
        si{1} = panel_islands_gen(nps,r,shift,nc);
        si{1}.n = 16*sum(nps);
        si{1}.np = sum(nps);

        shift = [0.75,-0.25]; r = 0.1;
        si{2} = panel_islands_gen(nps,r,shift,nc);
        si{2}.n = 16*sum(nps);
        si{2}.np = sum(nps);

        shift = [2.25,0.25]; r = 0.1;
        si{3} = panel_islands_gen(nps,r,shift,nc);
        si{3}.n = 16*sum(nps);
        si{3}.np = sum(nps);

        shift = [4,0.75]; r = 0.1;
        si{4} = panel_islands_gen(nps,r,shift,nc);
        si{4}.n = 16*sum(nps);
        si{4}.np = sum(nps);
    case 2
        nps = k*npsi;   % number of panels along each side
        si = {};
        shift = [-0.75,-0.75]; r = 0.1;
        si{1} = panel_islands_gen(nps,r,shift,nc);
        si{1}.n = 16*sum(nps);
        si{1}.np = sum(nps);

        shift = [0.75,-0.25]; r = 0.1;
        si{2} = panel_islands_gen(nps,r,shift,nc);
        si{2}.n = 16*sum(nps);
        si{2}.np = sum(nps);
        
        nps = k*[4;5;4;4;3;4;3;4;3;4;4;4];   % number of panels along each side
        r = 0.08;
        shift = [0.2,0.1];
        si{3} = panel_tank_gen(nps,r,shift,nc);    % number of panels at corner   
        si{3}.n = 16*sum(nps);
        si{3}.np = sum(nps);
        
end
        
