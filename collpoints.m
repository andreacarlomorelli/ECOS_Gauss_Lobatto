function [b, b_u] = collpoints(np, h, x, x_old, j, auxdata)
%
% collpoints
% 
% FUNCTION DESCRIPTION: It computes the equality dynamics constraints at
% the Gauss - Lobatto collocation points.
% 
% INPUTS:
% 
% PHI:     Matrix that finds the state at the collocation points  
%
% PHI_u:   Matrix that finds the control at the collocation points
%
% x:       State vector
%
% x_old:   State vector at the previous iteration
%
% PHI_p:   Derivative of the matrix that finds the state at the
%          collocation points                                      
%
% j:       Current segment of the Gauss - Lobatto method
%
% OUTPUTS:
%
% b        Vector of the constant terms for the state
%
% b_old    Old iteration's vector of the constant terms for the state
% 
% b_u      Vector of the constant terms for the control
%
% Delta:    Matrix of constraints
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 27/11/2020
%

% Sates and controls
n = 7;
m = 4;

% Useful parameters
c = auxdata.engine.c;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

Bm = B(c, ve, V0);

% Vector of constants terms: first part, 1 : np
b_old = zeros(2*np, n);
b(1 : np , :) = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);
b_old(1 : np,:) = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);

% Vector of constant terms: second part, np + 1 : 2 * np
x_n = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);
x_n_old = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);

u_n = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, n + 1 : n + m);
u_n_old = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, n + 1 : n + m);

for i = 1 : np
    b(np + i, :) = 0.5 * h * (f(x_n_old(i,:)) + A(x_n_old(i,:)) * ...
        (x_n(i,:)' - x_n_old(i,:)') + Bm*u_n(i,:)');
    b_old(np + i, :) = 0.5 * h * f_quad(x_n_old(i,:), u_n_old(i,:), c, ve, V0);
end

% Vector of constant terms: control
b_u(1 : np,:) = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, n + 1 : n + m);

end
