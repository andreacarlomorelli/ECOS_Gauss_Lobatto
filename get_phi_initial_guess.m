function [paraECOS, paraGL, paraSCP, auxdata] = get_phi_initial_guess(paraECOS, paraGL, paraSCP, auxdata, x)

% Calling useful parameters
n = 7;
m = 4;

ni = auxdata.ni;
nf = auxdata.nf;

% ECOS parameters
objective_c = paraECOS.objective_c;

% Gauss-Lobatto parameters
nc = paraGL.nc;
np = paraGL.np;
Ni = paraGL.Ni;
h = paraSCP.h;

% SCP parameters
muc = paraSCP.muc;
lambda = paraSCP.lambda;

c = auxdata.engine.c;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

% Upper and lower bounds on state and [taux, tauy, tauw]
G_upperb_bounds = paraECOS.G_upperb_bounds;
G_lowerb_bounds = paraECOS.G_lowerb_bounds;

h_upperb_bounds = paraECOS.h_upperb_bounds;
h_lowerb_bounds = paraECOS.h_lowerb_bounds;

% PHI matrices
PHI_p_matrix = sparse(auxdata.phi.PHI_p_matrix);
PHIu_matrix = sparse(auxdata.phi.PHIu_matrix);
PHIn_matrix_1 = sparse(auxdata.phi.PHIn_matrix_1);
PHIpn_matrix_1 = sparse(auxdata.phi.PHIpn_matrix_1);
PHIun_matrix_1 = sparse(auxdata.phi.PHIun_matrix_1);
PHIn_matrix_end = sparse(auxdata.phi.PHIn_matrix_end);
PHIpn_matrix_end = sparse(auxdata.phi.PHIpn_matrix_end);
PHIun_matrix_end = sparse(auxdata.phi.PHIun_matrix_end);

% Baseline  PHI matrices
PHI_c = sparse(paraGL.PHI_c); 
PHI_u = sparse(paraGL.PHI_u);

% Solution components' length
len_vect = paraECOS.len_vect;
x_len = len_vect(1);
u_len = len_vect(2);
virtual_tau_len = len_vect(5);
virtual_ctrl_len = len_vect(6);
sol_len = len_vect(11);

% HC_V is the vector of EQUALITY CONSTRAINTS of the original nonconvex,
% nonlinear problem
hc_v = zeros(1, virtual_ctrl_len + virtual_tau_len + ni + nf);

% Construction of the vector b_assembly
b_assembly = zeros(2*np*n*Ni,1);
bu_assembly = zeros(np*m*Ni,1);
f_quad_assembly = zeros(n*nc*Ni,1);

for j = 1 : Ni
    b_vect = zeros(2*np,n);
    x_n = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, 1 : n);
    u_n = x((j - 1) * (np - 1) + 1 : j * (np - 1) + 1, n + 1 : n + m);
    b_vect(1 : np,:) = x_n;
    bu_vect = u_n;
    for i = 1 : np
        b_vect(np + i, :) = 0.5 * h * (f_quad(x_n(i,:), u_n(i,:), c, ve, V0));
    end
    b_assembly((j-1)*2*np*n + 1 : ...
        (j-1)*2*np*n + 2*np*n) = reshape(b_vect,2*np*n,1);
    bu_assembly((j-1)*np*m + 1 : ...
        (j-1)*np*m + np*m) = reshape(bu_vect,np*m,1);
    for i = 1 : nc
        f_quad_assembly((j-1)*nc*n + (i-1)*n + 1 : (j-1)*nc*n + i*n) = ...
            f_quad(PHI_c(i,:)*b_vect, PHI_u(i,:)*bu_vect, c, ve, V0);
    end
end

J0 = objective_c(1 :  x_len + u_len)'*reshape(x, x_len + u_len, 1);

% Boundary collocation points: nonlinear dynamics equality constraints
hc_v(1 : n) = PHIpn_matrix_1*b_assembly(1 : 2*np*n) - ...
    h * 0.5 * f_quad(PHIn_matrix_1*b_assembly(1 : 2*np*n), PHIun_matrix_1*bu_assembly(1 : np*m), c, ve, V0);
hc_v(n+1 : 2*n) = PHIpn_matrix_end*b_assembly(end - 2*np*n + 1 : end) - ...
    h * 0.5 * f_quad(PHIn_matrix_end*b_assembly(end - 2*np*n + 1 : end), PHIun_matrix_end*bu_assembly(end - np*m + 1 : end), c, ve, V0);

% Intermidiate collocation points: nonlinear dynamics equality constraints
hc_v(2*n + 1 : virtual_ctrl_len) = PHI_p_matrix*b_assembly - ...
    h * 0.5 * f_quad_assembly;

% Boundary collocation points: nonlinear equality thrust constraints
tau_x_unit = [1 0 0 0];
tau_x_cell = repmat({tau_x_unit},1,Ni*nc);
tau_x_select = blkdiag(tau_x_cell{:});
tau_y_unit = [0 1 0 0];
tau_y_cell = repmat({tau_y_unit},1,Ni*nc);
tau_y_select = blkdiag(tau_y_cell{:});
tau_w_unit = [0 0 1 0];
tau_w_cell = repmat({tau_w_unit},1,Ni*nc);
tau_w_select = blkdiag(tau_w_cell{:});
tau_unit = [0 0 0 1];
tau_cell = repmat({tau_unit},1,Ni*nc);
tau_select = blkdiag(tau_cell{:});

hc_v(virtual_ctrl_len + 1 : virtual_ctrl_len + 1) = ...
    x(1, n+1)^2 + x(1, n+2)^2 + x(1, n+3)^2 - x(1, n+m)^2;
hc_v(virtual_ctrl_len + 2 : virtual_ctrl_len + 2) = ...
    x(end, n+1)^2 + x(end, n+2)^2 + x(end, n+3)^2 - x(end, n+m)^2;

% Intermidiate collocation points: nonlinear equality thrust constraints
hc_v(virtual_ctrl_len + 3 : virtual_ctrl_len + virtual_tau_len) = ...
    (tau_x_select*PHIu_matrix*bu_assembly).^2 + ...
    (tau_y_select*PHIu_matrix*bu_assembly).^2 + ...
    (tau_w_select*PHIu_matrix*bu_assembly).^2 - ...
    (tau_select*PHIu_matrix*bu_assembly).^2;

% Evaluation of the boundary conditions conditions constraints
hc_v(end-(ni+nf)+1 : end) = [x(1,1 : n) - paraSCP.x0(paraSCP.e,:), x(end,1 : n-1) - ...
    auxdata.xf];

% GC_V is the vector of INEQUALITY CONSTRAINTS of the original
% nonconvex, nonlinear problem
gc_v = zeros(1, 2*x_len + 2*u_len);

x_vect = zeros(sol_len,1);
x_vect(1 : x_len + u_len) = reshape(x, x_len + u_len, 1);

gc_v(1 : x_len + u_len) = G_upperb_bounds*x_vect - h_upperb_bounds;
gc_v(x_len + u_len + 1 : 2*x_len + 2*u_len) = ...
    G_lowerb_bounds*x_vect - h_lowerb_bounds;

paraECOS.phi_obj_0 = J0 + muc*norm(hc_v,1) + lambda*sum(max(0,gc_v));
paraECOS.phi_obj_hat_0 = paraECOS.phi_obj_0;

end

