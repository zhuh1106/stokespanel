function [ss,N] = polygon2(x,t,N,r)
% smooth polygon generating function
%
% x:   dim N x 2, N is the number of points to define a polygon,
%      need to be in order       
% r:   radius


[polySize, ~] = size(x);
% r = 0.1;

qtype = 'p';
qntype = 'C';
bdType = 'polygon';
idNum = 'middle';
% N = 4*48*polySize;

f = @(s) 0;
    i = 0;
    x1 = x(mod(i-1,polySize)+1,:);
    x2 = x(mod(i,polySize)+1,:);
    x3 = x(mod(i+1,polySize)+1,:);
    x4 = x(mod(i+2,polySize)+1,:);
%     tt = linspace(i*(2*pi)/polySize,(i+1)*(2*pi)/polySize,4);
    tt = [t(3*i+1); t(3*i+2); t(3*i+3); t(3*i+4)]';
    Zt = cInfBoundary2([x1; x2; x3; x4], r, tt, bdType, idNum);
    f = @(s) f(s) + (s<=tt(end)).*Zt(s);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);


for i = 1:polySize - 2
    x1 = x(mod(i-1,polySize)+1,:);
    x2 = x(mod(i,polySize)+1,:);
    x3 = x(mod(i+1,polySize)+1,:);
    x4 = x(mod(i+2,polySize)+1,:);
%     tt = linspace(i*(2*pi)/polySize,(i+1)*(2*pi)/polySize,4);
    tt = [t(3*i+1); t(3*i+2); t(3*i+3); t(3*i+4)]';
    Zt = cInfBoundary2([x1; x2; x3; x4], r, tt, bdType, idNum);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);
    f = @(s) f(s) + (s>tt(1)).*(s<=tt(end)).*Zt(s);
end

    i = polySize - 1;
    x1 = x(mod(i-1,polySize)+1,:);
    x2 = x(mod(i,polySize)+1,:);
    x3 = x(mod(i+1,polySize)+1,:);
    x4 = x(mod(i+2,polySize)+1,:);
%     tt = linspace(i*(2*pi)/polySize,(i+1)*(2*pi)/polySize,4);
    tt = [t(3*i+1); t(3*i+2); t(3*i+3); t(3*i+4)]';
    Zt = cInfBoundary2([x1; x2; x3; x4], r, tt, bdType, idNum);
    f = @(s) f(s) + (s>tt(1)).*Zt(s);


%% assembly
ss.Z = @(s) f(s);
% ss.Z = Z2;

%% quadrature pts get  
[ss, N, np] = quadr_pan(ss, N, qtype, qntype);

%% figure of cInf approx to triangle
% figure()
% % plot(real(ss.x),imag(ss.x),'.')
% plot(real(ss.x),imag(ss.x))
% daspect([1 1 1])
% 
% figure()
% plot(linspace(0,2*pi,N),ss.cur)
