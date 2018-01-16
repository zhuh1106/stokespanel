function z=It(x,p)
%
% It interpolates a periodic function at the points p. Multiple functions
% can be interpolated at a time by making x a matrix where each column
% corresponds to a periodic function. All functions must be on the interval 
% (0,2*pi] where only the right endpoint is included. To generate an 
% interpolation matrix, input the identity matrix for x.
%
% Example:
%   t=linspace(0,2*pi,65)'; t=t(2:end);
%   x=[exp(cos(t)),exp(sin(t))]; % functions to interpolate
%   p=rand(12,1); % points where the functions will be interpolated
%   z=It(x,p);

if nargin==0&&nargout==0
    STR=dbstack;
    STR=input(['Do you want to test ',STR.name,'.m? Use y/n: '],'s');
    if STR=='y'
        testcode
        return
    end
end

m=size(x,1);
X=fft(x([end,1:end-1],:));
K=[0:ceil(m/2)-1,-floor(m/2):-1];
z=real(exp(1i*p*K)*X/m);
end

function testcode
t=linspace(0,2*pi,65)';
t=t(2:end);
x=[exp(cos(t)),exp(sin(t))]; % functions to interpolate
p=rand(12,1); % points where the functions will be interpolated
EYE=eye(64); % identity matrix
z1=It(x,p);
error1=max(max(abs(z1-[exp(cos(p)),exp(sin(p))]))); % compare with known values
z2=It(EYE,p)*x;
error2=max(max(abs(z1-z2)));


% Test 1: Comparing with analytic solution
if error1<1e-14
    disp('Test 1: Success')
else
    disp('Test 1: Failure')
end

% Test 2: Comparing z1 and z2
if error2<1e-14
    disp('Test 2: Success')
else
    disp('Test 2: Failure')
end
end