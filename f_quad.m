function f = f_quad(x, u, c, ve, V0)
%
% f_quad
% 
% FUNCTION DESCRIPTION: Nonlinear, nonconvex dynamics vector
% 
% INPUTS:
% 
% x:    Vector of the states    [1x7]
%
% u:    Vector of the control variables    [1x4]
%
% OUTPUTS:
%
% f:    Vector of the dynamics  [7x1]
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

r = sqrt(x(1).^2 + x(2).^2 + x(3).^2);

f = [x(4);...
    x(5);...
    x(6);...
    -1./r.^3.*x(1) + c*u(1);...
    -1./r.^3.*x(2) + c*u(2);...
    -1./r.^3.*x(3) + c*u(3);...
    -c/(ve/V0)*u(4)];

end