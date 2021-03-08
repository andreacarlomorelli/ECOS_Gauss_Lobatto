function [paraECOS, paraGL, paraSCP, auxdata] = get_constant_matrices(paraECOS, paraGL, paraSCP, auxdata)
%
% get_constant_matrices
% 
% FUNCTION DESCRIPTION: this function computes the constant parts of the
% matrices G, A and of the vectors h, b
% 
% INPUTS/OUTPUTS:
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
N = paraGL.N;

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
virtual_tau_len = len_vect(5);
virtual_ctrl_len = len_vect(6);
aux_virtual_ctrl_len = len_vect(7);
aux_virtual_tau_len = len_vect(8);
slack_ineq_len = len_vect(9);
aux_trust_x_len = len_vect(10);
sol_len = len_vect(11);

% Call PHI matrices and vectors
PHIpn_matrix_1 = auxdata.phi.PHIpn_matrix_1;
PHIun_matrix_1 = auxdata.phi.PHIun_matrix_1;
PHI_p_matrix = auxdata.phi.PHI_p_matrix;
PHIu_matrix = auxdata.phi.PHIu_matrix;
PHIpn_matrix_end = auxdata.phi.PHIpn_matrix_end;
PHIun_matrix_end = auxdata.phi.PHIun_matrix_end;
weights_c_matrix = auxdata.trans.weights_c_matrix;
weights_n_matrix = auxdata.trans.weights_n_matrix;

% Call T and Tu
T_cost = auxdata.trans.T_cost;
Tu = auxdata.trans.Tu;

%%%%%% 1 : Adding slack variables to the A matrix %%%%%%

% Matrices initialization
A_dynamics_1 = zeros(n, sol_len);
A_dynamics_c = zeros(x_len_c, sol_len);
A_dynamics_end = zeros(n, sol_len);

% Matrix definition
A_dynamics_1(1 : n, u_len+x_len+1 : u_len+x_len+n*np) = eye(n,n*np);
A_dynamics_end(1 : n, u_len+x_len+1 : u_len+x_len+n*np) = eye(n,n*np);
for i = 1 : Ni
    A_dynamics_c((i-1)*n*nc + 1 : i*n*nc, u_len+x_len+(i-1)*n*np + 1 : u_len+x_len+i*n*np) = eye(n*nc,n*np);
end

% Building constant parts of the A matrix
A_dynamics_1(:, :) = A_dynamics_1(:, :) + PHIpn_matrix_1*T_cost(1 : 2*np*n, :);
A_dynamics_c(:,:) = A_dynamics_c + PHI_p_matrix*T_cost;
A_dynamics_end(:, :) = A_dynamics_end(1 : n, :) + PHIpn_matrix_end*T_cost(end - 2*np*n + 1 : end, :);

% Saving A matrix components
paraECOS.A_dynamics_1_cost = A_dynamics_1;
paraECOS.A_dynamics_c_cost = A_dynamics_c;
paraECOS.A_dynamics_end_cost = A_dynamics_end;

%%%%%% 2 : Relaxed boundary constraints of x0 and xf %%%%%%
paraECOS.G_upperb_x0 = zeros(n, sol_len);
paraECOS.G_lowerb_x0 = zeros(n, sol_len);
paraECOS.G_upperb_xf = zeros(n-2, sol_len);
paraECOS.G_lowerb_xf = zeros(n-2, sol_len);

% Relaxed initial/final position
x0 = paraSCP.x0(paraSCP.e,:);
xf = paraSCP.xf(paraSCP.e,:);

paraECOS.h_upperb_x0 = x0' + auxdata.eps_x0;
paraECOS.h_lowerb_x0 = - x0' + auxdata.eps_x0;
paraECOS.h_upperb_xf = xf' + auxdata.eps_xf;
paraECOS.h_lowerb_xf = - xf' + auxdata.eps_xf;

% Matrices for initial/final position
paraECOS.G_upperb_x0(1:n, 1:n) = eye(n);
paraECOS.G_lowerb_x0(1:n, 1:n) = -1*eye(n);
paraECOS.G_upperb_xf(1:n-2, n*(N-1)+1:n*N-2) = eye(n-2);
paraECOS.G_lowerb_xf(1:n-2, n*(N-1)+1:n*N-2) = -1*eye(n-2);

%%%%%% 3 : Bounds on x and u %%%%%%
G_upperb_bounds = zeros(x_len+u_len, sol_len);
G_upperb_bounds(1:end, 1:x_len+u_len) = eye(x_len+u_len);
G_lowerb_bounds = zeros(x_len+u_len, sol_len);
G_lowerb_bounds(1:end, 1:x_len+u_len) = -1*eye(x_len+u_len);

paraECOS.G_upperb_bounds = G_upperb_bounds;
paraECOS.G_lowerb_bounds = G_lowerb_bounds;

% Auxiliary variables to define lower and upper bounds
upperb_h_x_tmp = auxdata.bounds.ub(1 : n);
lowerb_h_x_tmp = -1*auxdata.bounds.lb(1 : n);

upperb_h_x = repmat(upperb_h_x_tmp,1,N);
lowerb_h_x = repmat(lowerb_h_x_tmp,1,N);
upperb_h_x = upperb_h_x(:);
lowerb_h_x = lowerb_h_x(:);

upperb_h_u_tmp = auxdata.bounds.ub(n + 1 : n + m);
lowerb_h_u_tmp = -1*auxdata.bounds.lb(n + 1 : n + m);
upperb_h_u = repmat(upperb_h_u_tmp,1,N);
lowerb_h_u = repmat(lowerb_h_u_tmp,1,N);
upperb_h_u = upperb_h_u(:);
lowerb_h_u = lowerb_h_u(:);

% Combine bounds
h_upperb_bounds = [upperb_h_x; upperb_h_u];
h_lowerb_bounds = [lowerb_h_x; lowerb_h_u];

paraECOS.h_upperb_bounds = h_upperb_bounds;
paraECOS.h_lowerb_bounds = h_lowerb_bounds;

%%%%%% 4 : Constraint for tau magnitude %%%%%%
G_upperb_tau = zeros(Ni*nc + 2, sol_len);
G_lowerb_tau = zeros(Ni*nc + 2, sol_len);

% Selection of z
z_unit = [0 0 0 0 0 0 1];
z_cell = repmat({z_unit},1,Ni*nc);
paraECOS.z_matrix_select = blkdiag(z_cell{:});
paraECOS.z_unit = z_unit;

% Selection of tau
tau_unit = [0 0 0 1];
tau_cell = repmat({tau_unit},1,Ni*nc);
tau_matrix_select = blkdiag(tau_cell{:});

% Adding constant part of the constraints on tau

G_upperb_tau(1, :) = tau_unit*PHIun_matrix_1*Tu(1 : np*m, :);
G_upperb_tau(2:end-1,:) = tau_matrix_select*PHIu_matrix*Tu;
G_upperb_tau(end, :) = tau_unit*PHIun_matrix_end*Tu(end - np*m + 1 : end, :);

G_lowerb_tau(1,:) =  -tau_unit*PHIun_matrix_1*Tu(1 : np*m, :);
G_lowerb_tau(2:end-1,:) = -tau_matrix_select*PHIu_matrix*Tu;
G_lowerb_tau(end,:) = -tau_unit*PHIun_matrix_end*Tu(end - np*m + 1 : end, :);

paraECOS.G_upperb_tau_cost = G_upperb_tau;
paraECOS.G_lowerb_tau_cost = G_lowerb_tau;

%%%%%% 5 : Constraint for virtual control %%%%%%

% Index of the solution/optimization vector where aux_vc starts
aux_vc_start = x_len + u_len + virtual_ctrl_len;

G_vc_slack1 = zeros(aux_virtual_ctrl_len, sol_len);
G_vc_slack2 = zeros(aux_virtual_ctrl_len, sol_len);

% First additional constraint vc - aux_vc <= 0
G_vc_slack1(:, x_len+u_len+1:x_len+u_len+virtual_ctrl_len) = eye(virtual_ctrl_len); % Virtual_ctrl
G_vc_slack1(:, aux_vc_start+1:aux_vc_start+aux_virtual_ctrl_len) = -1*eye(aux_virtual_ctrl_len); % aux_virtual_ctrl

paraECOS.G_vc_slack1 = G_vc_slack1;

% Second additional constraint -vc - aux_vc <= 0
G_vc_slack2(:, x_len+u_len+1:x_len+u_len+virtual_ctrl_len) = -1*eye(virtual_ctrl_len); % virtual_ctrl
G_vc_slack2(:, aux_vc_start+1:aux_vc_start+aux_virtual_ctrl_len) = -1*eye(aux_virtual_ctrl_len); % aux_virtual_ctrl

paraECOS.G_vc_slack2 = G_vc_slack2;

% Corresponding h vectors
paraECOS.h_vc_slack1 = zeros(aux_virtual_ctrl_len,1);
paraECOS.h_vc_slack2 = zeros(aux_virtual_ctrl_len,1);

%%%%%% 6 : Constraint for virtual tau %%%%%%

% Index of the solution/optimization vector where aux_tau starts
aux_tau_start = x_len + u_len + virtual_ctrl_len + aux_virtual_ctrl_len + virtual_tau_len;

G_tau_slack1 = zeros(aux_virtual_tau_len, sol_len);
G_tau_slack2 = zeros(aux_virtual_tau_len, sol_len);

% First additional constraint virtual_tau - aux_tau <= 0
G_tau_slack1(:, x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+1:x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+virtual_tau_len) = eye(virtual_tau_len); % virtual_tau
G_tau_slack1(:, aux_tau_start+1:aux_tau_start+aux_virtual_tau_len) = -1*eye(aux_virtual_tau_len);
paraECOS.G_tau_slack1 = G_tau_slack1;

% Second additional constraint - virtual_tau - aux_tau <= 0
G_tau_slack2(:, x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+1:x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+virtual_tau_len) = -1*eye(virtual_tau_len); % virtual_tau
G_tau_slack2(:, aux_tau_start+1:aux_tau_start+aux_virtual_tau_len) = -1*eye(aux_virtual_tau_len);

paraECOS.G_tau_slack2 = G_tau_slack2;

% Corresponding h vectors
paraECOS.h_tau_slack1 = zeros(aux_virtual_tau_len,1);
paraECOS.h_tau_slack2 = zeros(aux_virtual_tau_len,1);

% Nonnegativity of slack_ineq (= bounds of this variable)
slack_ineq_start = x_len + u_len;

G_slack_ineq_bounds = zeros(slack_ineq_len, sol_len);

G_slack_ineq_bounds(:, slack_ineq_start+1:slack_ineq_start+slack_ineq_len) = -1*eye(slack_ineq_len);

paraECOS.G_slack_ineq_bounds = G_slack_ineq_bounds;

paraECOS.h_slack_ineq_bounds = zeros(slack_ineq_len, 1);

%%%%%% 7 : Second Order Cone Constraint %%%%%%
socc_matrix = zeros(u_len_c + 2*m,sol_len);

submatrix = [0 0 0 -1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0];   % Original sub matrix that is to be repeated on the diagonal
socc = repmat({submatrix},1,Ni*nc);
socc_matrix(1:m,:) = submatrix*PHIun_matrix_1*Tu(1 : np*m, :);
socc_matrix(m+1:end-m,:) = blkdiag(socc{:})*PHIu_matrix*Tu;
socc_matrix(end-m+1:end,:) = submatrix*PHIun_matrix_end*Tu(end - np*m + 1 : end, :);
paraECOS.G_socc_control = socc_matrix;

paraECOS.h_socc_control = zeros(u_len_c + 2*m, 1);

%%%%%% 8 : Constant parts of the trust region %%%%%%

% matrix/vector for the constraint 1^T * aux_trust <= rk
paraECOS.G_trust = zeros(N, sol_len);
paraECOS.h_trust = zeros(N,1);

% matrices/vectors for -aux_trust-x_k <= -x_k_ref and x_k - aux_trust <= x_k_ref
paraECOS.G_trust_x_slack1 = zeros(aux_trust_x_len, sol_len);
paraECOS.G_trust_x_slack2 = zeros(aux_trust_x_len, sol_len);
paraECOS.h_trust_x_slack1 = zeros(aux_trust_x_len,1);
paraECOS.h_trust_x_slack2 = zeros(aux_trust_x_len,1);

aux_trust_x_start = x_len + u_len + virtual_ctrl_len + ...
    aux_virtual_ctrl_len + virtual_tau_len + aux_virtual_tau_len;

% Trust region (inequality constraints)
% || xk - xk_ref ||_1 <= rk

for i=1:N % we consider all points

    % -aux_trust-x_k <= -x_k_ref
    paraECOS.G_trust_x_slack1(n*(i-1)+1: n*i, n*(i-1)+1: n*i) = -1*eye(n); % x
    paraECOS.G_trust_x_slack1(n*(i-1)+1: n*i, aux_trust_x_start+n*(i-1)+1: aux_trust_x_start+n*i) = -1*eye(n); % aux_x_trust
    
    % x_k - aux_trust <= x_k_ref
    paraECOS.G_trust_x_slack2(n*(i-1)+1: n*i, n*(i-1)+1: n*i) = 1*eye(n); % x
    paraECOS.G_trust_x_slack2(n*(i-1)+1: n*i, aux_trust_x_start+n*(i-1)+1: aux_trust_x_start+n*i) = -1*eye(n); % aux_x_trust
    
    % 1^T * aux_trust <= rk
    paraECOS.G_trust(i, aux_trust_x_start+n*(i-1)+1:aux_trust_x_start+n*i) = ones(1,n);
    
end

% Objective function c_objective
objective_c = zeros(sol_len,1);

aux_slack_ineq_start = x_len + u_len;
objective_c(aux_slack_ineq_start+1:aux_slack_ineq_start+slack_ineq_len) = paraSCP.muc * ones(slack_ineq_len, 1);

% Virtual control
aux_vc_start = x_len + u_len + virtual_ctrl_len;
objective_c(aux_vc_start+1:aux_vc_start+aux_virtual_ctrl_len) = paraSCP.muc * ones(aux_virtual_ctrl_len, 1);

% Virtual tau
aux_tau_start = x_len + u_len + virtual_ctrl_len + aux_virtual_ctrl_len;

objective_c(aux_tau_start+1:aux_tau_start+aux_virtual_tau_len) = paraSCP.lambda * ones(aux_virtual_tau_len, 1);

% Maximize final mass
obj_coll_points = sum(weights_c_matrix*PHIu_matrix*Tu,1);
obj_coll_points = obj_coll_points(x_len+1:x_len+u_len)';
obj_node_points = sum(weights_n_matrix*Tu,1);
obj_node_points = obj_node_points(x_len+1:x_len+u_len)';

objective_c(x_len+1:x_len+u_len) = h * 0.5 * (obj_coll_points + obj_node_points);

paraECOS.objective_c = objective_c;


end

