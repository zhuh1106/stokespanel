function st = panel_tank_gen(nps,r,shift,nc)

addpath(genpath('.')) % adds subdirectories to path
warning off % turns off warning from gmres

% generate tank geometry
x1 = [0,-0.2]; x2 = [0.75,-1]; x3 = [3.25,-1]; x4 = [4,-0.2]; x5 = [2.5,-0.2];
x6 = [2.5,0.25]; x7 = [4,0.25]; x8 = [4,0.55]; x9 = [2.5,0.55]; x10 = [2.5,1];
x11 = [1.5,1]; x12 = [1.5,-0.2];
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];
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
[st,~] = polygon_test2(x,r,N,t);