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

f = [x(:,1)*x(:,3)/x(:,4);...
    x(:,1)*x(:,5)/x(:,4);...
    x(:,4) - 1/(x(:,1)*x(:,4));...
    - x(:,3);...
    - x(:,2)/(x(:,1)^2*x(:,4));...
    x(:, 1)/x(:,4);...
    0];

end

