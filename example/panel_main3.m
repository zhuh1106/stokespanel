function panel_main3()

close all
clear all

xLim = 6;
r = 0.1;
nc = 2;
k = 2;
nps = k*[4;5;5;8;5;5;4;5;4;5;4;5;5;8;5;5];   % number of panels along each side
%% get geometry
sc = panel_confinedgeo_gen2(xLim,nps,r,nc);
sc.n = numel(sc.x);
sc.np = sum(nps);
N = numel(sc.x);

%% set up rhs (boundary condition) and solve for density
f1 = [(1-imag(sc.x).^2).*(abs(imag(sc.x))<1).*(real(sc.x)<(-xLim+1));zeros(N,1)];
f2 = [(1-(2*(abs(imag(sc.x))-1)).^2).*(abs(imag(sc.x))>0.5).*(abs(imag(sc.x))<1.5).*...
    (real(sc.x)>(xLim-1));zeros(N,1)];
f = f1 + f2;
figure(),quiver(real(sc.x),imag(sc.x),f(1:end/2),f((end/2+1):end)); hold on; 
plot(real(sc.x),imag(sc.x),'.'); axis equal 
idx1 = find(f1>0); idx2 = find(f2>0);
sum(sc.w(idx1).*f1(idx1)), sum(sc.w(idx2).*f2(idx2)),
%% set up islands
option = 3;
si = islands(k,nc,option);

numOfI = 40;


%% set up target point
nx = 160; gx = ((1:nx)/nx*2-1)*xLim; ny = 80; gy = ((1:ny)/ny*2-1)*3;
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
x = [-xLim xLim]; y = [-3,3];
u0 = 1; u1c = min(max(real(u),-u0),u0); u2c = min(max(imag(u),-u0),u0);
starty = gy; startx = gx(1)*ones(size(starty));
imagesc(x,y,sqrt((u1c.*ii).^2+(u2c.*ii).^2)); colormap(jet(256));hold on
plot(real(sc.x),imag(sc.x),'k'); hold on
quiver(xx,yy, u1c.*ii,u2c.*ii, 0.8); hold on
streamline(xx,yy,u1c.*ii,u2c.*ii,startx,starty);


% quiver(xx,yy, real(u).*ii,imag(u).*ii, 0.8);
axis equal
c = [0.75 0.75 0.75];
for i = 1:numOfI
    plot(real(si{i}.x),imag(si{i}.x),'k');
    fill(real(si{i}.x),imag(si{i}.x),c);
end
fillout(real(sc.x),imag(sc.x),[-xLim-1 xLim+1 -5 5],1.2*c);

function si = islands(k,nc,option)

switch option
    case 1
        nps = k*[8;3;8;3];   % number of panels along each side
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
        nps = k*[8;3;8;3];   % number of panels along each side
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
    case 3
        scale = 1/2; theta = pi/6;
        nps = k*[3;3;3];   % number of panels along each side
        si = {};
        Shift = [-7/24,5/8; -7/24,5/8-3/4; -7/24,5/8-3/4*2; -7/24,5/8-3/4*3;...
                -7/24,-4-5/8+3/4*3; -7/24,-4-5/8+3/4*2; -7/24,-4-5/8+3/4; -7/24,-4-5/8];
        r = 0.1;
        for j = 1:8
            for i = 1:8
                shift = Shift(i,:);
                shift(1) = shift(1) + 3/4*(j-1);
                si{i+(j-1)*8} = panel_islands_gen2(nps,r,shift,nc,scale,theta);
                si{i+(j-1)*8}.n = 16*sum(nps);
                si{i+(j-1)*8}.np = sum(nps);
            end
        end   
        
end
        
