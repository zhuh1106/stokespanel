function [xt, yt] = cInfFun(theta, r)
% c infinity function to smooth a corner

alpha = 12 * cot(theta) / exp(2);   % derive stretching factor to match angle

xt = @(t) alpha * t * r / ( 1/2 + alpha^2 * exp(2) / 24);   
                                    % parametric form for x coordinate
                                    % multiplied by similarity factor
                                                 
yt = @(t) (mypsi(t) ./ (mypsi(t) + mypsi(1-t)) ...
            - (1/2-alpha^2*exp(2)/24) ) * r / ( 1/2 + alpha^2 * exp(2) / 24);
                                    % paramatric form for y coordinate
                                    % multiplied by similarify factor
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
