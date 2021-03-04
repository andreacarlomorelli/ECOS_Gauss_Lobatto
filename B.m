function [B] = B(c, ve, V0)
%
% B
% 
% FUNCTION DESCRIPTION: this function finds the B matrix of the linearized
% dynamics of the problem x_dot = f(x_star) + A(x_star)*(x - x_star) + B*u.
% 
% INPUTS:
% 
% c         Problem parameter
% 
% ve        Exhaust velocity
%
% V0        Velocities adimensionalization factor
%
% OUTPUTS:
%
% B         B matrix
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 26/02/2021
%

% B matrix
B = [0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    c 0 0 0;...
    0 c 0 0;...
    0 0 c 0;...
    0 0 0 -c/(ve/V0)];

end

