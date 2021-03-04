function [dx] = dyn_rec(t, x, c, ve, V0, tvect, tau, taux, tauy, tauw)
%
% dyn_rec
% 
% FUNCTION DESCRIPTION: 
% 
% INPUTS:
% 
% t:    Time instant                                    [1x1]
%
% x:    State vector                                    [7x1]
%
% tvect:  Vector of nondimensional time                 [1 x variable]
%
% tau:  Control variable tau                            [1 x variable]
%
% taux: Control variable taux                           [1 x variable]
%
% tauy: Control variable tauy                           [1 x variable]
%
% tauw: Control variable tauw                           [1 x variable]

% OUTPUTS:
%
% dx:   Derivatives of the state
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

tau = interp1(tvect,tau,t,'spline');
taux = interp1(tvect,taux,t,'spline');
tauy = interp1(tvect,tauy,t,'spline');
tauw = interp1(tvect,tauw,t,'spline');

r = sqrt(x(1,:)^2 + x(2,:)^2 + x(3,:)^2);

dx = zeros(7,1);

dx(1,:) = x(4,:);

dx(2,:) = x(5,:);

dx(3,:) = x(6,:);

dx(4,:) = - 1/r^3*x(1,:) + c*taux;

dx(5,:) = - 1/r^3*x(2,:) + c*tauy;

dx(6,:) = - 1/r^3*x(3,:) + c*tauw;

dx(7,:) = -c/(ve/V0)*tau;

end
