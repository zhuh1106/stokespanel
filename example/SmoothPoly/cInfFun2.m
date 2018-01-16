function [xt, yt] = cInfFun2(theta, r)
% a different choice of c infinity function

alpha = pi*cot(theta);   % derive stretching factor to match angle

x = @(t) alpha * t * r / ( alpha / 2 * sqrt(1+(alpha/pi)^2));  
% x = @(t) alpha * (t - pi/2) * r / ( alpha * pi / 2 * sqrt(1+alpha^2));   
                                    % parametric form for x coordinate
                                    % multiplied by similarity factor
                                                 
y = @(t) ( cos(pi*t) + alpha^2 / (2*pi) )...
            * r / ( alpha / 2 * sqrt(1+(alpha/pi)^2));
                                    % paramatric form for y coordinate
                                    % multiplied by similarify factor
xt = @(t) -cos(theta)*x(1/2 - t) + sin(theta)*y(1/2 - t);
yt = @(t) sin(theta)*x(1/2 - t) + cos(theta)*y(1/2 - t);


                                    

% function test_cInfFun                              
% theta = pi/6;
% r = 0.5;   
% [xt, yt] = cInfFun( theta, r);
% ss.Z = @(t) xt(t/(4*pi)) + 1i*yt(t/(4*pi));
% 
% [ss, N, np] = quadr_pan(ss, N, qtype, qntype);
% 
% figure()
% plot(real(ss.x),imag(ss.x),'*')
% xlim([0 0.6])
% ylim([0 0.6])
