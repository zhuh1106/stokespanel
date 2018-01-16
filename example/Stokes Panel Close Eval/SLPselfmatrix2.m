function S=SLPselfmatrix2(s)
%
% This function generates the single-layer self interaction matrix using
% the Hybrid Gauss-Trapezoidal rule. The function only works on one vesicle
% at a time.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   t=t(2:end); % Note that only the right endpoint is included.
%   s.x=[cos(t)+sin(2*t)/4,sin(t)/2]*[1;1i];
%   S=SLPselfmatrix(s);

if nargin==0&&nargout==0
    STR=dbstack;
    STR=input(['Do you want to test ',STR.name,'.m? Use y/n: '],'s');
    if STR=='y'
        testcode
        return
    end
end

ss=[real(s.x([end,1:end])),imag(s.x([end,1:end]))];
N=size(ss,1); t=linspace(0,2*pi,N)'; C=tril(ones(N-1)); C=logical([C;~C]);
Y1=ss(2:end,1)*ones(1,N-1); Y1=[Y1;Y1]; Y1=reshape(Y1(C),[N-1,N-1]).';
Y2=ss(2:end,2)*ones(1,N-1); Y2=[Y2;Y2]; Y2=reshape(Y2(C),[N-1,N-1]).';
Y1=[Y1,Y1(:,1)]; Y1=[Y1(end,:);Y1]; Y2=[Y2,Y2(:,1)]; Y2=[Y2(end,:);Y2];
Y=[Y1,Y2]; Y=[Y,Dtp(Y)]; G=@(x,y)SLPG(y,ss);

% Alpert quadrature is done here.
S=Alpert(G,t,Y,1);

S11=S(1:end/3,:); S11(:,end)=S11(:,1)+S11(:,end); S11=S11(2:end,2:end);
S12=S(end/3+1:2*end/3,:); S12(:,end)=S12(:,1)+S12(:,end);
S12=S12(2:end,2:end); S22=S(2*end/3+1:end,:);
S22(:,end)=S22(:,1)+S22(:,end); S22=S22(2:end,2:end); C=C(:,end:-1:1);
S11=[S11,S11].'; S11=reshape(S11(C),[N-1,N-1]).';
S12=[S12,S12].'; S12=reshape(S12(C),[N-1,N-1]).';
S22=[S22,S22].'; S22=reshape(S22(C),[N-1,N-1]).';
S=[S11,S12;S12,S22];
end

function z=SLPG(y,s)
%
% The single layer kernel is defined here.
%

N=size(y,1); Y1=y(:,1:end/4); Y2=y(:,end/4+1:end/2);
DY1=y(:,end/2+1:3*end/4); DY2=y(:,3*end/4+1:end);
R1=ones(N,1)*s(:,1)'-Y1; R2=ones(N,1)*s(:,2)'-Y2;
RHO2=R1.^2+R2.^2; LOGRHO=-log(sqrt(RHO2)); SA=sqrt(DY1.^2+DY2.^2);
z11=1/4/pi*(LOGRHO+(R1.*R1)./RHO2).*SA;
z12=1/4/pi*((R1.*R2)./RHO2).*SA;
z22=1/4/pi*(LOGRHO+(R2.*R2)./RHO2).*SA;
z=[z11,z12,z22];
end

function z=Alpert(G,t,y,d)

m=length(t); if nargin==4&&d==1, L0=t(end)-t(1); L=L0/2; m0=m; m=(m+1)/2; 
else L=t(end)-t(1); end
% Alpert nodes and weights for 8th order integration.
v=[6.531815708567918e-3,9.086744584657729e-2,3.967966533375878e-1,...
    1.027856640525646e0,1.945288592909266e0,2.980147933889640e0,...
    3.998861349951123e0]';
u=[2.462194198995203e-2,1.701315866854178e-1,4.609256358650077e-1,...
    7.947291148621894e-1,1.008710414337933e0,1.036093649726216e0,...
    1.004787656533285e0]';
x=[2.087647422032129e-1,9.786087373714483e-1,1.989541386579751e0,3e0]';
w=[5.207988277246498e-1,9.535038018555888e-1,1.024871626402471e0,...
    1.000825744017291e0]';
g=11;n=round(m-9);h=L/(n+8);a=5;

if nargin==4&&d==1
    IT=Itp(eye(m0),[v*h;L-x*h]);
    tg=[v*h;L-x*h;t(a+1:a+n)];tg=[tg;L0-tg];
    if isempty(y), Geval=h*G(tg).'; else yud=y(end:-1:1,:);
    yg=[IT*y;y(a+1:a+n,:);IT*yud;yud(a+1:a+n,:)]; Geval=h*G(tg,yg).'; end
    W=(ones(size(Geval,1),1)*[u;w]');
    zL=(W.*Geval(:,1:g))*IT;
    zR=(W.*Geval(:,end/2+1:end/2+g))*IT;
    zL(:,a+1:a+n)=zL(:,a+1:a+n)+Geval(:,g+1:end/2);
    zR(:,a+1:a+n)=zR(:,a+1:a+n)+Geval(:,end/2+g+1:end);
    z=zL+zR(:,end:-1:1);
else
    IT=Itp(eye(m),[v*h;L-x*h]);
    tg=[v*h;L-x*h;t(a+1:a+n)]; if isempty(y), Geval=h*G(tg).'; else
    yg=[IT*y;y(a+1:a+n,:)]; Geval=h*G(tg,yg).'; end
    z=((ones(size(Geval,1),1)*[u;w]').*Geval(:,1:g))*IT;
    z(:,a+1:a+n)=z(:,a+1:a+n)+Geval(:,g+1:end);
end
end

function Dx=Dtp(x,d)
if nargin==1
    d=1;
end
x=x(2:end,:);
Dx=Dt(x,d);
Dx=Dx([end,1:end],:);
end

function xout=Itp(x,p)
x=x(2:end,:);
xout=It(x,p);
end

function testcode
mu=0.7; %viscosity
a=linspace(0,2*pi,257)';
a=a(2:end);
s.x=[cos(a) sin(a)+cos(2*a)/2]*[1;1i]; % vesicle position
s=quadr(s); % Computes properties of the vesicle
tau=[exp(sin(a)),exp(cos(a))]*[1;1i]; % density
TAU=[real(tau);imag(tau)];% density as a real vector

% Test 2: Comparing self interaction with StokesScloseeval
z1=SLPselfmatrix(s)*TAU/mu;
z1=z1(1:end/2)+1i*z1(end/2+1:end); % include tau/2 for PV limit
z2=StokesScloseeval(s.x+1e-14*s.nx,s,tau,'e')/mu; % close evaluation scheme
error1=max(abs(z1-z2));
if error1<1e-11
    disp('Test 1: Success')
else
    disp('Test 1: Failure')
end
end