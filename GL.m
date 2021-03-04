function [paraGL] = GL(paraGL)
%
% GL
% 
% FUNCTION DESCRIPTION: this function computes the points, the weights and
% the matrices related to a given arbitrary order Gauss - Lobatto method.
% 
% INPUTS:
% 
% N:    Order of the Gauss - Lobatto method   [1x1]
%
% OUTPUTS:
%
% para:    Structure containing different outputs
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

ng = paraGL.ng;
[x, w, ~]=lglnodes(ng - 1);

% Gauss-Lobatto parameters
paraGL.nc = (ng - 1)*0.5; % Number of collocation points per segment
paraGL.np = ng - paraGL.nc; % Number of nodes per segment
paraGL.N = paraGL.np + (paraGL.np - 1)*(paraGL.Ni - 1); % Total number of nodes

% Definition of the nodes and collocation points
paraGL.tau_g = x(end : -2 : 1);
paraGL.zeta = x(end - 1 : -2 : 2);

% Definition of the weights at nodes and collocation points
paraGL.weights_n = w(end : -2 : 1);
paraGL.weights_c = w(end - 1 : -2 : 2);

% Matrices definition
paraGL.A_gl = A_gauss(ng, paraGL.tau_g);
paraGL.Z_gl = Z_matr_gauss(ng,paraGL.zeta);
paraGL.A_gl_u = A_gauss_u(ng, paraGL.tau_g);
paraGL.Zdot_gl = Zdot_matr_gauss(ng,paraGL.zeta);
paraGL.Z_gl_u = Z_matr_gauss_u(ng,paraGL.zeta);
paraGL.PHI_c = paraGL.Z_gl*pinv(paraGL.A_gl);
paraGL.PHI_pc = paraGL.Zdot_gl*pinv(paraGL.A_gl);
paraGL.Z_gl_n = Z_matr_gauss_n(ng,paraGL.tau_g);
paraGL.Zdot_gl_n = Zdot_matr_gauss_n(ng,paraGL.tau_g);
paraGL.PHI_n = paraGL.Z_gl_n*pinv(paraGL.A_gl);
paraGL.PHI_pn = paraGL.Zdot_gl_n*pinv(paraGL.A_gl);
paraGL.PHI_u = paraGL.Z_gl_u*pinv(paraGL.A_gl_u);

% Auxiliary functions definition
function [x,w,P]=lglnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde 
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods. 
%
% Reference on LGL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
             
end

w=2./(N*N1*P(:,N1).^2);
end
function A = A_gauss(ng, tau)
%
% A_gauss
% 
% FUNCTION DESCRIPTION: It computes the matrix of coefficients of the
% polynomial that interpolates the state at the nodal points. 
% The state is interpolated with a polynomial of order ng + 1 where ng 
% is the order of the Gauss - Lobatto method.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method                 [1x1]
%
% tau:   Position of the nodes inside the interval (-1,1)    [1 x variable]
%
% OUTPUTS:
%
% A:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

A = zeros( ng + 1, ng + 1);
A(1 : ( ng + 1 )*0.5, 1) = 1;
A(( ng + 1 )*0.5 + 1 : end, 1) = 0;
A(( ng + 1 )*0.5 + 1 : end, 2) = 1;

for i = 1 : size(A,1)*0.5
    for j = 2 : size(A,2)
        A(i,j) = tau(i)^(j-1);
    end
end

for i = size(A,1)*0.5 + 1 : size(A,1)
    for j = 3 : size(A,2)
        A(i,j) = (j-1)*tau(i - size(A,1)*0.5)^(j-2);
    end
end

end
function A = A_gauss_u(ng, tau)
%
% A_gauss_u
% 
% FUNCTION DESCRIPTION: It computes the matrix of coefficients of the
% polynomial that interpolates the control at the nodal points. The control is interpolated
% with a polynomial of order np where np is the number of nodes per segment.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method               [1x1]
%
% tau:   Position of the nodes inside the interval (-1,1)  [1 x variable]
%
% OUTPUTS:
%
% A:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

A = zeros((ng + 1)*0.5, (ng + 1)*0.5);
A(:, 1) = 1;

for i = 1 : size(A,1)
    for j = 2 : size(A,2)
        A(i,j) = tau(i)^(j-1);
    end
end

end
function Z = Z_matr_gauss(ng, zeta)
%
% Z_matr_gauss
% 
% FUNCTION DESCRIPTION: It computes the matrix of coefficients of the
% polynomial that interpolates the state at the collocation points. 
% The state is interpolated with a polynomial of order ng + 1 where ng 
% is the order of the Gauss - Lobatto method.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method            [1x1]
%
% zeta:   Position of the collocation points            [1 x variable]
%         inside the interval (-1,1)    
%
% OUTPUTS:
%
% Z:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Z = zeros((ng - 1)*0.5, ng + 1);
Z(:,1) = 1;

for i = 1 : size(Z,1)
    for j = 2 : size(Z,2) 
        Z(i,j) = zeta(i)^(j-1);
    end
end

end
function Z = Z_matr_gauss_n(ng, tau)
%
% Z_matr_gauss_n
% 
% FUNCTION DESCRIPTION: It computes the matrix of coefficients of the
% polynomial that interpolates the state at the initial and final nodes. 
% The state is interpolated with a polynomial of order ng + 1 where ng 
% is the order of the Gauss - Lobatto method.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method            [1x1]
%
% tau:   Position of the nodal points                   [1 x variable]
%        inside the interval (-1,1)    
%
% OUTPUTS:
%
% Z:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Z = zeros((ng + 1)*0.5, ng + 1);
Z(:,1) = 1;

for i = 1 : size(Z,1)
    for j = 2 : size(Z,2) 
        Z(i,j) = tau(i)^(j-1);
    end
end

end
function Z = Z_matr_gauss_u(ng, zeta)
%
% Z_matr_gauss_u
% 
% FUNCTION DESCRIPTION: It computes the matrix of coefficients of the
% polynomial that interpolates the control at the collocation points. 
% The state is interpolated with a polynomial of order np where np 
% is the number of nodes per segment.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method            [1x1]
%
% tau:   Position of the nodal points                   [1 x variable]
%        inside the interval (-1,1)    
%
% OUTPUTS:
%
% Z:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Z = zeros((ng - 1)*0.5, (ng + 1)*0.5);
Z(:,1) = 1;

for i = 1 : size(Z,1)
    for j = 2 : size(Z,2) 
        Z(i,j) = zeta(i)^(j-1);
    end
end

end
function Zdot = Zdot_matr_gauss(ng, zeta)
%
% Zdot_matr_gauss
% 
% FUNCTION DESCRIPTION: It computes the derivative of the matrix of coefficients of the
% polynomial that interpolates the state at the collocation points. 
% The state is interpolated with a polynomial of order ng + 1 where ng 
% is the order of the Gauss - Lobatto method.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method            [1x1]
%
% zeta:   Position of the collocation points            [1 x variable]
%        inside the interval (-1,1)    
%
% OUTPUTS:
%
% Z:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Zdot = zeros((ng - 1)*0.5, ng + 1);
Zdot(:,1) = 0;
Zdot(:,2) = 1;

for i = 1 : size(Zdot,1)
    for j = 3 : size(Zdot,2) 
        Zdot(i,j) = (j-1)*zeta(i)^(j-2);
    end
end

end
function Zdot = Zdot_matr_gauss_n(ng, zeta)
%
% Zdot_matr_gauss_n
% 
% FUNCTION DESCRIPTION: It computes the derivative of the matrix of coefficients of the
% polynomial that interpolates the state at the initial and final nodes. 
% The state is interpolated with a polynomial of order ng + 1 where ng 
% is the order of the Gauss - Lobatto method.
% 
% INPUTS:
% 
% ng:    Order of the Gauss - Lobatto method            [1x1]
%
% zeta:   Position of the collocation points            [1 x variable]
%        inside the interval (-1,1)    
%
% OUTPUTS:
%
% Z:    Matrix of coefficients
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 19/09/2020
%

Zdot = zeros((ng + 1)*0.5, ng + 1);
Zdot(:,1) = 0;
Zdot(:,2) = 1;

for i = 1 : size(Zdot,1)
    for j = 3 : size(Zdot,2) 
        Zdot(i,j) = (j-1)*zeta(i)^(j-2);
    end
end

end

end

