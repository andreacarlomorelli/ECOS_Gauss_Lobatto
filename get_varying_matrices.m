function [paraGL, paraSCP, paraECOS, auxdata] = get_varying_matrices(x_old, paraGL, paraSCP, paraECOS, auxdata)
%
% get_varying_matrices
% 
% FUNCTION DESCRIPTION: this function computes the varying parts of the
% matrices G, A and of the vectors b, h.
% 
% INPUTS:
%
% x_old          Reference solution                          
% 
% paraECOS:      Structure with the ECOS parameters 
%
% paraGL:        Structure with the GL parameters 
%
% paraSCP:       Structure with the SCP parameters 
%
% auxdata:       Structure with the auxiliary parameters
%
% OUTPUTS:                          
% 
% paraECOS:      Structure with the ECOS parameters 
%
% paraGL:        Structure with the GL parameters 
%
% paraSCP:       Structure with the SCP parameters 
%
% auxdata:       Structure with the auxiliary parameters
% 
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 25/02/2021
%

% Sates and controls
n = 7;
m = 4;

% ECOS parameters
len_vect = paraECOS.len_vect;

% GL parameters
nc = paraGL.nc;
np = paraGL.np;
Ni = paraGL.Ni;

% SCP parameters
h = paraSCP.h;

% Useful parameters
c = auxdata.engine.c;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

% Length of the variables
x_len = len_vect(1);
u_len = len_vect(2);
x_len_c = len_vect(3);
u_len_c = len_vect(4);
sol_len = len_vect(11);

% PHI matrices
PHI_c = sparse(paraGL.PHI_c);
PHI_matrix = sparse(auxdata.phi.PHI_matrix);
PHI_p_matrix = sparse(auxdata.phi.PHI_p_matrix);
PHIu_matrix = sparse(auxdata.phi.PHIu_matrix);
PHIpn_matrix_1 = sparse(auxdata.phi.PHIpn_matrix_1);
PHIn_matrix_1 = sparse(auxdata.phi.PHIn_matrix_1);
PHIun_matrix_1 = sparse(auxdata.phi.PHIun_matrix_1);
PHIpn_matrix_end = sparse(auxdata.phi.PHIpn_matrix_end);
PHIn_matrix_end = sparse(auxdata.phi.PHIn_matrix_end);
PHIun_matrix_end = sparse(auxdata.phi.PHIun_matrix_end);

% Transformation matrices and vectors
T = sparse(auxdata.trans.T);
Tu = sparse(auxdata.trans.Tu);
T_cost = sparse(auxdata.trans.T_cost);
As = paraECOS.As;
Bs = paraECOS.Bs;

% Vectors b_old_assembly, b_cost_assembly and f_dyn and matrix A_dyn
A_dyn = zeros(x_len_c);
B_dyn = zeros(x_len_c, u_len_c);
f_dyn = zeros(x_len_c,1);
b_old_assembly = zeros(2*x_len_c,1);
b_cost_assembly = zeros(2*x_len_c,1);

for j = 1 : Ni
    b_cost = zeros(2*np,n);
    b_old = zeros(2*np,n);
    x_n_old = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);
    u_n_old = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, n + 1 : n + m);
    b_old(1 : np,:) = x_old((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);
    b_cost(1 : np,:) = zeros(np,n);
    for i = 1 : np
        b_old(np + i, :) = 0.5 * h * f_quad(x_n_old(i,:), u_n_old(i,:), c, ve, V0);
        b_cost(np + i, :) = 0.5* h * (f(x_n_old(i,:)) - ...
            A(x_n_old(i,:)) * x_n_old(i,:)');
    end
    b_old_assembly((j-1)*2*np*n + 1 : ...
        (j-1)*2*np*n + 2*np*n) = reshape(b_old,2*np*n,1);
    b_cost_assembly((j-1)*2*np*n + 1 : ...
        (j-1)*2*np*n + 2*np*n) = reshape(b_cost,2*np*n,1);
    for i = 1 : nc
        A_par = A(PHI_c(i,:)*b_old);
        B_par = B(PHI_c(i,:)*b_old, c, ve, V0);
        f_par = f(PHI_c(i,:)*b_old);
        % A_dyn
        A_dyn((j-1)*n*nc + (i-1)*n + 1 : (j-1)*n*nc + (i-1)*n + n, ...
            (j-1)*n*nc + (i-1)*n + 1 : (j-1)*n*nc + (i-1)*n + n) = A_par;
        % B_dyn
        B_dyn((j-1)*n*nc + (i-1)*n + 1 : (j-1)*n*nc + (i-1)*n + n, ...
            (j-1)*m*nc + (i-1)*m + 1 : (j-1)*m*nc + (i-1)*m + m) = B_par;
        % f_lyn
        f_dyn((j-1)*n*nc + (i-1)*n + 1 : (j-1)*n*nc + (i-1)*n + n) = f_par;  
    end
end

A_dyn = sparse(A_dyn);
B_dyn = sparse(B_dyn);

paraECOS.b_cost_assembly = b_cost_assembly;
paraECOS.b_old_assembly = b_old_assembly;

% Assembled A matrix and vector b
A_dynamics_1_cost = sparse(paraECOS.A_dynamics_1_cost);
A_dynamics_c_cost = sparse(paraECOS.A_dynamics_c_cost);
A_dynamics_end_cost = sparse(paraECOS.A_dynamics_end_cost);

% Selecting only the varying part of T
T_var = T - T_cost;

A_dynamics_1 = A_dynamics_1_cost + PHIpn_matrix_1*T_var(1 : 2*np*n, :) - h * 0.5 * sparse(As(:,:,1))*PHIn_matrix_1*T(1 : 2*np*n, :) - h * 0.5 * (sparse(Bs(:,:,1))*PHIun_matrix_1*Tu(1 : np*m, :));
A_dynamics_c = A_dynamics_c_cost + PHI_p_matrix*T_var - h * 0.5 * (A_dyn*PHI_matrix*T) - h * 0.5 * (B_dyn*PHIu_matrix*Tu);
A_dynamics_end = A_dynamics_end_cost + PHIpn_matrix_end*T_var(end - 2*np*n + 1 : end, :) - h * 0.5 * sparse(As(:,:,end))*PHIn_matrix_end*T(end - 2*np*n + 1 : end, :) - h * 0.5 * (sparse(Bs(:,:,end))*PHIun_matrix_end*Tu(end - np*m + 1 : end, :));

paraECOS.A_dynamics = [A_dynamics_1; A_dynamics_c; A_dynamics_end];

paraECOS.A_dyn_c = A_dynamics_c;

b_dynamics_1 = -PHIpn_matrix_1*b_cost_assembly(1 : 2*np*n) + ...
    h * 0.5 * (f(x_old(1,1:n)) + A(x_old(1,1:n))*(PHIn_matrix_1*...
    b_cost_assembly(1 : 2*np*n) - x_old(1, 1 : n)'));

b_dynamics = -PHI_p_matrix*b_cost_assembly + h * 0.5 * (f_dyn + ...
    sparse(A_dyn)*PHI_matrix*(b_cost_assembly - b_old_assembly));

b_dynamics_end = -PHIpn_matrix_end*b_cost_assembly(end - 2*np*n + 1 : end) ...
    + h * 0.5 * (f(x_old(end,1:n)) + A(x_old(end,1:n))*(PHIn_matrix_end*...
    b_cost_assembly(end - 2*np*n + 1 : end) - x_old(end, 1 : n)'));

paraECOS.b_dynamics = [b_dynamics_1; b_dynamics; b_dynamics_end];

paraECOS.b_dyn_c = b_dynamics;

% Constraint for tau magnitude
G_upperb_tau_cost = sparse(paraECOS.G_upperb_tau_cost);
G_lowerb_tau_cost = sparse(paraECOS.G_lowerb_tau_cost);

G_upperb_tau = zeros(Ni*nc + 2, sol_len);
G_lowerb_tau = zeros(Ni*nc + 2, sol_len);

% Selection of z
z_unit = paraECOS.z_unit;
z_matrix_select = sparse(paraECOS.z_matrix_select);

x_old_vect = zeros(sol_len, 1);
x_old_vect(1 : x_len) = reshape(x_old(:, 1 : n)',x_len,1);
x_old_vect(x_len+1 : x_len+u_len) = reshape(x_old(:,n+1:n+m)',u_len,1);
x_old_c = PHI_matrix*T*x_old_vect;
z_old_c = z_matrix_select*x_old_c(:);

% First and last collocation points
z_old_1 = x_old(1,n);
z_old_end = x_old(end,n);

D_1 = exp(-z_old_1)*(1 + z_old_1); 
D = exp(-z_old_c).*(1 + z_old_c);
D_end = exp(-z_old_end)*(1 + z_old_end); 

G_upperb_tau(1, :) = G_upperb_tau_cost(1, :) + exp(-z_old_1)*z_unit*PHIn_matrix_1*T(1: 2*np*n,:);
G_upperb_tau(2:end-1,:) = G_upperb_tau_cost(2:end-1,:) + exp(-z_old_c).*z_matrix_select*PHI_matrix*T;
G_upperb_tau(end, :) = G_upperb_tau_cost(end, :) + exp(-z_old_end)*z_unit*PHIn_matrix_end*T(end - 2*np*n + 1 : end, :);

G_lowerb_tau(1,:) = G_lowerb_tau_cost(1,:) - exp(-z_old_1).*z_unit*PHIn_matrix_1*T(1: 2*np*n,:);
G_lowerb_tau(2:end-1,:) = G_lowerb_tau_cost(2:end-1,:) - exp(-z_old_c).*z_matrix_select*PHI_matrix*T;
G_lowerb_tau(end,:) = G_lowerb_tau_cost(end,:) -  exp(-z_old_end).*z_unit*PHIn_matrix_end*T(end - 2*np*n + 1 : end, :);

paraECOS.G_upperb_tau = G_upperb_tau;
paraECOS.G_lowerb_tau = G_lowerb_tau;

% Associated constant vectors
paraECOS.h_upperb_tau = [D_1 + 1e-7 - z_unit*PHIn_matrix_1*b_cost_assembly(1 : 2*np*n); ...
    D - z_matrix_select*PHI_matrix*b_cost_assembly; ...
    D_end - z_unit*PHIn_matrix_end*b_cost_assembly(end - 2*np*n + 1 : end)];

paraECOS.h_lowerb_tau = [1 - z_unit*PHIn_matrix_1*b_cost_assembly(1 : 2*np*n); ...
    1 - z_matrix_select*PHI_matrix*b_cost_assembly; ...
   1 - z_unit*PHIn_matrix_end*b_cost_assembly(end - 2*np*n + 1 : end)];

end

