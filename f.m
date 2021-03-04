function f = f(x)
%
% f
% 
% FUNCTION DESCRIPTION: this function finds the vector f of the linearized
% dynamics of the problem x_dot = f(x_star) + A(x_star)*(x - x_star) + B*u.
% 
% INPUTS:
% 
% x:    Vector of the states    [variable x 7]
%
% OUTPUTS:
%
% f:    Vector of the dynamics  [7x1]
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

r = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

f = [x(:,4);...
    x(:,5);...
    x(:,6);...
    -1./r.^3.*x(:,1);...
    -1./r.^3.*x(:,2);...
    -1./r.^3.*x(:,3);...
    0];

end

