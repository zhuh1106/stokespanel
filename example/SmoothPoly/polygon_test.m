function [ss,N] = polygon_test(x,r,k)
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
N = k*48*polySize;

f = @(s) 0;
for i = 0:polySize - 1
    x1 = x(mod(i,polySize)+1,:);
    x2 = x(mod(i+1,polySize)+1,:);
    x3 = x(mod(i+2,polySize)+1,:);
    x4 = x(mod(i+3,polySize)+1,:);
    t = linspace(i*(2*pi)/polySize,(i+1)*(2*pi)/polySize,4);
    Zt = cInfBoundary2([x1; x2; x3; x4], r, t, bdType, idNum);
%     ss.Z = @(s) (s<=pi).*Zt1(s) + (s>pi).*Zt2(s);
    f = @(s) f(s) + (s>t(1)).*(s<=t(end)).*Zt(s);
end


%% assembly
ss.Z = @(s) (s==0).*f(2*pi)+(s>0).*f(s);
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
