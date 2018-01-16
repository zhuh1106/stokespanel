function si = panel_islands_gen(nps,r,shift,nc)

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

%% set up source

% generate microfludic chip geometry
x1 = [-2+shift,2]; x2 = [-2+shift,-2]; x3 = [-1+shift,-2]; x4 = [-1+shift,2];
x = [x1;x2;x3;x4];
np = sum(nps); theta = 2*pi/np;
[polySize, ~] = size(x);
t = zeros(polySize,4);
t(1,:) = [0,nc*theta,nps(1)/np*2*pi-nc*theta,nps(1)/np*2*pi];
for i = 2:polySize
    t(i,:) = [t(i-1,end),t(i-1,end)+nc*theta,t(i-1,end)+(nps(i)-nc)*theta,t(i-1,end)+nps(i)*theta];
end

N = 16*np;
[si,~] = polygon_test2(x,r,N,t);