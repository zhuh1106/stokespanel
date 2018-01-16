% triangle need to be consistent with raidus r
function [result,R] = radius_test(x1,x2,x3,r)
% close all
% clear all
% x1 = [0,0];
% x2 = [1,1];
% x3 = [1,0];

result = 0;
vec = [x1;x2;x3];
area = 1/2*norm(cross(vec(:,1),vec(:,2)));
R = 2*area/(norm(x2-x1)+norm(x3-x1)+norm(x3-x2));
if R > r
    result = 1;
end

% N = 128;
% 
% t = linspace(0,2*pi,10);
% 
% ss = triangle_test(x1, x2, x3, t, N, r);