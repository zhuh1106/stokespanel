close all
clear all

% randomly generate 4 points
x = randi([-10 10],1,4);

% coordinates of 3 corners
x1 = [0,0];
x2 = [x(1),x(2)];
x3 = [x(3),x(4)];

% number of quadrature points
N = 9*16*2;     % 3 side has at least 9 panels, each panel has 16 pts, 
                % multiplied by any number (here it is 2)

% starting and ending points in parametrization for this piecewise defined smooth triangle
t = linspace(0,2*pi,10);

% pre-assigned radius
r = 0.2;
[result,R] = radius_test(x1,x2,x3,r);

if result == 1
    N = 16*9*6;
    
    [ss,N] = triangle_test(x1, x2, x3, t, N, r);
    figure()
    plot(linspace(0,2*pi,N),ss.cur)
    title('curvature')
    hold on
    plot(t,zeros(size(t)),'*')
    hold off
    
    fprintf('smooth triangle generated\n')
else
    fprintf('radius is too large\n')
    
    [ss,N] = triangle_test(x1, x2, x3, t, N, 0.9*R);
    figure()
    plot(linspace(0,2*pi,N),ss.cur)
    title('curvature')
    hold on
    plot(t,zeros(size(t)),'*')
    hold off
end