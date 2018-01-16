close all
clear all
r = 0.2;

x1 = [-100,1]; x2 = [-100,-1]; x3 = [100,-1]; x4 = [100,1];
x = [x1;x2;x3;x4]; k = 1;

[ss,N] = polygon_test(x,r,k);

x = [0,3;1,1;3,2;2,-1;3,-3;0,-2;-3,-3;-2,-1;-3,2;-1,1];
[ss,N] = polygon_test(x,r,1);
ss;