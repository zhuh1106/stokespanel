function [su, N] = pipe_test(xu, r, k)
% smooth pipe generating function
%
% xu:  dim N x 2, N is the number of points to define one side of a pipe,
%      need to be in order       
% r:   radius

[uWallSize, ~] = size(xu); uWallSize = uWallSize - 1;

qtype = 'p';
qntype = 'C';
bdType = 'polygon';
N = k*48*uWallSize;

%% predefine function, and then add each piece
f = @(s) 0;

%% first piece of wall
idNum = 'first';
x1 = xu(1,:); x2 = xu(2,:); x3 = xu(3,:);
% t = linspace(0,2*pi/uWallSize,3);
t = [0,4*pi/(3*uWallSize),2*pi/uWallSize];
Zt = cInfBoundary2([x1; x2; x3], r, t, bdType, idNum);
f = @(s) f(s) + (s>=t(1)).*(s<t(end)).*Zt(s);

%% middle pieces of wall
idNum = 'middle';
for i = 2:uWallSize - 1 % 2nd to 2nd to last
    x1 = xu(i-1,:); x2 = xu(i,:); x3 = xu(i+1,:); x4 = xu(i+2,:);
    t = linspace((i-1)*(2*pi)/uWallSize,i*(2*pi)/uWallSize,4);
    Zt = cInfBoundary2([x1; x2; x3; x4], r, t, bdType, idNum);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);
    f = @(s) f(s) + (s>=t(1)).*(s<t(end)).*Zt(s);
end

%% last piece of wall
idNum = 'last';
x1 = xu(end-2,:); x2 = xu(end-1,:); x3 = xu(end,:);
% t = linspace((uWallSize-1)*(2*pi)/uWallSize,2*pi,3);
t = [(uWallSize-1)*(2*pi)/uWallSize,2*pi-4*pi/(3*uWallSize),2*pi];
Zt = cInfBoundary2([x1; x2; x3], r, t, bdType, idNum);
f = @(s) f(s) + (s>=t(1)).*(s<=t(end)).*Zt(s);


%% assembly
su.Z = @(s) f(s);
% ss.Z = Z2;

%% quadrature pts get  
[su, N, np] = quadr_pan(su, N, qtype, qntype);

%% figure of cInf approx to triangle
% figure()
% % plot(real(ss.x),imag(ss.x),'.')
% plot(real(su.x),imag(su.x),'.')
% daspect([1 1 1])
% 
% figure()
% plot(linspace(0,2*pi,N),su.cur)

