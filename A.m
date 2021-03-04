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

r = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

A = [0 0 0 1 0 0 0;...
    0 0 0 0 1 0 0;...
    0 0 0 0 0 1 0;...
   -(-2*x(:,1).^2 + x(:,2).^2 + x(:,3).^2)./r.^5, 3*x(:,1).*x(:,2)./r.^5, 3*x(:,1).*x(:,3)./r.^5 0 0 0 0;...
    3*x(:,1).*x(:,2)./r.^5, -(x(:,1).^2 - 2*x(:,2).^2 + x(:,3).^2)./r.^5, 3*x(:,2).*x(:,3)./r.^5 0 0 0 0;...
    3*x(:,1).*x(:,3)./r.^5, 3*x(:,2).*x(:,3)./r.^5, -(x(:,1).^2 + x(:,2).^2 - 2*x(:,3).^2)./r.^5 0 0 0 0;...
    0 0 0 0 0 0 0];
   
    

end

