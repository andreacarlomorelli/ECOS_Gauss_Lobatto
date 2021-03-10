function [dx] = dyn_rec(t, x, c, ve, V0, tvect, tau, taur, tauth, tauw)
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
taur = interp1(tvect,taur,t,'spline');
tauth = interp1(tvect,tauth,t,'spline');
tauw = interp1(tvect,tauw,t,'spline');

s = sqrt(x(1).^2 + x(2).^2);

dx = zeros(7,1);

dx(1,:) = x(1,:)*x(3,:)/x(4,:);

dx(2,:) = x(1,:)*x(5,:)/x(4,:);

dx(3,:) = - x(1,:)^2/(s^3*x(4,:)) + x(4,:) + c*x(1,:)/x(4,:)*taur;

dx(4,:) = -x(3) + c*x(1,:)/x(4,:)*tauth;

dx(5,:) = -x(1)*x(2)/(s^3*x(4)) + c*x(1,:)/x(4,:)*tauw;

dx(6,:) = x(1,:)/x(4,:);

dx(7,:) = -c/(ve/V0)*x(1,:)/x(4,:)*tau;

end
