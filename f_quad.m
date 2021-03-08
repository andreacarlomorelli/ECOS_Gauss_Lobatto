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

s = sqrt(x(1).^2 + x(2).^2);

f = [x(1)*x(3)/x(4);...
    x(1)*x(5)/x(4);...
    -x(1)^2/(s^3*x(4)) + x(4) + c*x(1)./x(4).*u(1);...
    -x(3) + c*x(1)./x(4).*u(2);...
    -x(1)*x(2)/(s^3*x(4)) + c*x(1)./x(4).*u(3);...
    x(1)./x(4);...
    - c/(ve/V0)*x(1)./x(4).*u(4)];

end