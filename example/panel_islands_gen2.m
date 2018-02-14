function si = panel_islands_gen2(nps,r,shift,nc,scale,theta)

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

if nargin < 6
    theta = 0;
end
if nargin < 5
    scale = 1;
end

%% set up source

% generate microfludic chip geometry
x1 = [-3,2]; x2 = [-2,1.5]; x3 = [-2,2.5];
xc = 1/3*(x1 + x2 + x3);
rmat = [cos(theta),-sin(theta);sin(theta),cos(theta)];
x = scale*(rmat*(([x1;x2;x3] - xc)'))'+xc;
shiftm = zeros(size(x));
shiftm(:,1) = shift(1);
shiftm(:,2) = shift(2);
x = x + shiftm;
np = sum(nps); theta = 2*pi/np;
[polySize, ~] = size(x);
t = zeros(polySize,4);
t(1,:) = [0,nc*theta,nps(1)/np*2*pi-nc*theta,nps(1)/np*2*pi];
for i = 2:polySize
    t(i,:) = [t(i-1,end),t(i-1,end)+nc*theta,t(i-1,end)+(nps(i)-nc)*theta,t(i-1,end)+nps(i)*theta];
end

N = 16*np;
[si,~] = polygon_test2(x,r,N,t);