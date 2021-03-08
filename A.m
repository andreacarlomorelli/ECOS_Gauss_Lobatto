function A = A(x)
%
% A
% 
% FUNCTION DESCRIPTION: this function finds the A matrix of the linearized
% dynamics of the problem x_dot = f(x_star) + A(x_star)*(x - x_star) + B*u.
% 
% INPUTS:
% 
% x:    Vector of the states    [variable x 7]
%
% OUTPUTS:
%
% A:    Jacobian of the states  [7x7]
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

s = sqrt(x(1)^2 + x(2)^2);

A = [x(3)/x(4) 0 x(1)/x(4) -x(1)*x(3)/x(4)^2 0 0 0;...
    x(5)/x(4) 0 0 -x(1)*x(5)/x(4)^2 x(1)/x(4) 0 0;...
    (x(1)^3 - 2*x(1)*x(2)^2)/(x(:,4)*s^5) 3*x(1)^2*x(2)/(x(:,4)*s^5) 0 1+x(1)^2/(s^3*x(:,4)^2) 0 0 0;...
    0 0 -1 0 0 0 0; ...
    -x(2)*(x(2)^2 - 2*x(1)^2)/(x(4)*s^5) -x(1)*(x(1)^2 - 2*x(2)^2)/(x(:,4)*s^5) 0 x(1)*x(2)/(x(4)^2*s^3) 0 0 0; ...
    1/x(4) 0 0 -x(1)/x(4)^2 0 0 0; ...
    0 0 0 0 0 0 0];
   
    

end

