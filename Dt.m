function z=Dt(x,d)
%
% This program computes the dth derivative of a periodic function on 
% the interval (0,2*pi]. If only one input is given, the function outputs
% the first derivative.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   t=t(2:end);
%   y=sin(t);
%   dy=Dt(y);
%   ddy=Dt(y,2);


if nargin==0&&nargout==0
    STR=dbstack;
    STR=input(['Do you want to test ',STR.name,'.m? Use y/n: '],'s');
    if STR=='y'
        testcode
        return
    end
end

[m,p]=size(x);
if nargin==1
    d=1;
end

X=fft(x([end,1:end-1],:));
K=1i*[0:ceil(m/2)-1,-floor(m/2):-1]'*ones(1,p);
Kd=K.^d;
if d==-1
   Kd(isinf(Kd))=0;
end
z=ifft(Kd.*X);
z=real(z([2:end,1],:));
if d==-1
    z=z-z(end);
    t=linspace(0,2*pi,m+1)';
    t=t(2:end);
    z=t*real(X(1,:))/m+z;
end
end

function testcode
n=128;
t=linspace(0,2*pi,n+1)';
t=t(2:end);
y=[exp(sin(t)),exp(cos(t))];
dy=Dt(y);
ddy=Dt(y,2);

% Test 1: Checking 1st derivative with analytic solution
error1=max(max(abs(dy-[exp(sin(t)).*cos(t),-exp(cos(t)).*sin(t)])));
if error1<1e-13
    disp('Test 1: Success')
else
    disp('Test 1: Failure')
end

% Test 2: Checking 2nd derivative with analytic solution
error2=max(max(abs(ddy-[exp(sin(t)).*(cos(t).^2-sin(t)),exp(cos(t)).*(sin(t).^2-cos(t))])));
if error2<1e-11
    disp('Test 2: Success')
else
    disp('Test 2: Failure')
end

% Test 3: Checking Integration (odd)
t=linspace(0,2*pi,64)';
t=t(2:end);
error3=max(abs(Dt(sin(t).^2,-1)-1/2*(t-sin(t).*cos(t))));
if error3<1e-14
    disp('Test 3: Success')
else
    disp('Test 3: Failure')
end

% Test 4: Checking Integration (even)
t=linspace(0,2*pi,65)';
t=t(2:end);
error4=max(abs(Dt(sin(t).^2,-1)-1/2*(t-sin(t).*cos(t))));
if error4<1e-14
    disp('Test 4: Success')
else
    disp('Test 4: Failure')
end

% Test 5: Checking Integration Vectorization
t=linspace(0,2*pi,65)';
t=t(2:end);
error5=max(max(abs(Dt([sin(t),cos(t)],-1)-[-cos(t),sin(t)])));
if error5<1e-14
    disp('Test 5: Success')
else
    disp('Test 5: Failure')
end
end