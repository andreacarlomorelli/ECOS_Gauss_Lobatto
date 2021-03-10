function [paraECOS, paraGL, paraSCP, x, x_old_out, time, J0] = ECOS_SCP(paraECOS, paraGL, paraSCP, auxdata)
%
% ECOS_SCP
% 
% FUNCTION DESCRIPTION: this function implements the Sequential Convex
% Programming Algorithm (SCP).
% 
% INPUTS:                          
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
% ecos_result:   Solution vector
%
% time:          Vector of times of each ECOS iteration
%
% AUTHOR: Andrea Carlo Morelli
%
% DATE: 25/02/2021
%

% Sates and controls
n = 7;
m = 4;

% Initial guess
x_old = paraSCP.x_old;

% ECOS parameters
len_vect = paraECOS.len_vect;

% GL parameters
nc = paraGL.nc;
N = paraGL.N;
Ni = paraGL.Ni;

% SCP parameters
h = paraSCP.h; lambda = paraSCP.lambda; muc = paraSCP.muc;
epsc = paraSCP.epsc; epsphi = paraSCP.epsphi; r_tr = paraSCP.r_tr;
rho0 = paraSCP.rho0; rho1 = paraSCP.rho1; rho2 = paraSCP.rho2;
alpha = paraSCP.alpha; alpha_max = paraSCP.alpha_max; beta = paraSCP.beta;
beta_max = paraSCP.beta_max; delta = paraSCP.delta; gamma = paraSCP.gamma;

% SCP initialization parameters
cmax = 1; dphi = 1; k = 0; phi_obj_old = paraECOS.phi_obj_0;
phi_hat_old = paraECOS.phi_obj_hat_0; ratio_muc = 1; ratio_lambda = 1;  

% Initializaing runtime vector
time = [];

% Initializing alpha and beta vectors
paraSCP.alpha_vect = [];
paraSCP.beta_vect = [];

% Length of the variables
x_len = len_vect(1);
u_len = len_vect(2);
virtual_tau_len = len_vect(5);
virtual_ctrl_len = len_vect(6);
aux_virtual_ctrl_len = len_vect(7);
aux_trust_x_len = len_vect(10);
sol_len = len_vect(11);

% Calling constant ECOS matrices and vectors

% Boundary Conditions
G_upperb_x0 = paraECOS.G_upperb_x0;
G_lowerb_x0 = paraECOS.G_lowerb_x0;
G_upperb_xf = paraECOS.G_upperb_xf;
G_lowerb_xf = paraECOS.G_lowerb_xf;

h_upperb_x0 = paraECOS.h_upperb_x0;
h_lowerb_x0 = paraECOS.h_lowerb_x0;
h_upperb_xf = paraECOS.h_upperb_xf;
h_lowerb_xf = paraECOS.h_lowerb_xf;

% Virtual control
G_vc_slack1 = paraECOS.G_vc_slack1;
G_vc_slack2 = paraECOS.G_vc_slack2;

h_vc_slack1 = paraECOS.h_vc_slack1;
h_vc_slack2 = paraECOS.h_vc_slack2;

% Virtual tau
G_tau_slack1 = paraECOS.G_tau_slack1;
G_tau_slack2 = paraECOS.G_tau_slack2;

h_tau_slack1 = paraECOS.h_tau_slack1;
h_tau_slack2 = paraECOS.h_tau_slack2;

% Second Order Cone Constraints
G_socc_control = paraECOS.G_socc_control;

h_socc_control = paraECOS.h_socc_control;

% Slack variables trust region
G_trust_x_slack1 = paraECOS.G_trust_x_slack1;
G_trust_x_slack2 = paraECOS.G_trust_x_slack2;

h_trust_x_slack1 = paraECOS.h_trust_x_slack1;
h_trust_x_slack2 = paraECOS.h_trust_x_slack2;

% Bounds of the slack variables
G_slack_ineq_bounds = paraECOS.G_slack_ineq_bounds;

h_slack_ineq_bounds = paraECOS.h_slack_ineq_bounds;

% Upper and lower bounds on state and [taux, tauy, tauw]
G_upperb_bounds = paraECOS.G_upperb_bounds;
G_lowerb_bounds = paraECOS.G_lowerb_bounds;

h_upperb_bounds = paraECOS.h_upperb_bounds;
h_lowerb_bounds = paraECOS.h_lowerb_bounds;

% Upper and lower bounds on tau
G_upperb_tau = paraECOS.G_upperb_tau;
G_lowerb_tau = paraECOS.G_lowerb_tau;

h_upperb_tau = paraECOS.h_upperb_tau;
h_lowerb_tau = paraECOS.h_lowerb_tau;

% Objective function
objective_c = paraECOS.objective_c;

% Trust region
G_trust = paraECOS.G_trust;

h_trust = paraECOS.h_trust;

% ECOS requires to give the sizes of the matrices in a struct
dimensions_cones.l = size(G_trust,1) + size(G_trust_x_slack1,1) + size(G_trust_x_slack1,1) + ...
    size(G_vc_slack1,1) + size(G_vc_slack2,1) + ...
    + size(G_tau_slack1,1) +  size(G_tau_slack2,1) + ...
    size(G_upperb_bounds,1) + size(G_lowerb_bounds,1) + ...
    size(G_upperb_tau,1) + size(G_lowerb_tau,1) + ...
    size(G_lowerb_xf,1) + size(G_upperb_xf,1) + ...
    size(G_lowerb_x0,1) + size(G_upperb_x0,1) + ...
    size(G_slack_ineq_bounds,1);

% SOOC (see also documentation):
q_total_tmp = size(G_socc_control,1) + 2*m;
q_per_constraint_tmp = q_total_tmp / (Ni*nc + 2);
q_control = q_per_constraint_tmp * ones(1,(Ni*nc + 2));
dimensions_cones.q = q_control;

% This part is only to optimize run time
h_trust_x_slack1_start = N;
h_trust_x_slack1_end = h_trust_x_slack1_start + aux_trust_x_len;
h_trust_x_slack2_start = h_trust_x_slack1_end;
h_trust_x_slack2_end = h_trust_x_slack2_start + aux_trust_x_len;

h_trust_indices = 1:N;
h_trust_x_slack1_indices = h_trust_x_slack1_start+1:h_trust_x_slack1_end;
h_trust_x_slack2_indices = h_trust_x_slack2_start+1:h_trust_x_slack2_end;

% ECOS opts
opts.VERBOSE = 0;

while cmax >= epsc || dphi >= epsphi
    % Get varying parts of matrices T and Tu
    [paraECOS, paraGL, paraSCP, auxdata] = get_varying_T_Tu(x_old, paraECOS, paraGL, paraSCP, auxdata);
    
    % Get varying parts of matrices A_ecos, G_ecos, h_ecos and b_ecos
    [paraGL, paraSCP, paraECOS, auxdata] = get_varying_matrices(x_old, paraGL, paraSCP, paraECOS, auxdata);
    
    % Calling varying ECOS matrices and vectors
    
    % Dynamics
    A_dynamics = paraECOS.A_dynamics;
    b_dynamics = paraECOS.b_dynamics;
    
    % Upper and lower bounds on state and [taux, tauy, tauw]
    G_upperb_bounds = paraECOS.G_upperb_bounds;
    G_lowerb_bounds = paraECOS.G_lowerb_bounds;
    
    h_upperb_bounds = paraECOS.h_upperb_bounds;
    h_lowerb_bounds = paraECOS.h_lowerb_bounds;
    
    % Upper and lower bounds on tau
    G_upperb_tau = paraECOS.G_upperb_tau;
    G_lowerb_tau = paraECOS.G_lowerb_tau;
    
    h_upperb_tau = paraECOS.h_upperb_tau;
    h_lowerb_tau = paraECOS.h_lowerb_tau;
    
    % Objective function
    objective_c = paraECOS.objective_c;
    
    % Trust region
    G_trust = paraECOS.G_trust;
    
    h_trust = paraECOS.h_trust;
    
    % Combine all matrices
    G_ecos = [sparse(G_trust); sparse(G_trust_x_slack1); sparse(G_trust_x_slack2); ...
        sparse(G_vc_slack1); sparse(G_vc_slack2); sparse(G_tau_slack1); sparse(G_tau_slack2); ...
        sparse(G_upperb_bounds); sparse(G_lowerb_bounds); sparse(G_lowerb_x0); sparse(G_upperb_x0);
        sparse(G_lowerb_xf); sparse(G_upperb_xf); sparse(G_lowerb_tau); ...
        sparse(G_upperb_tau); sparse(G_slack_ineq_bounds); sparse(G_socc_control)];
    
    % Combine vectors
    h_ecos = [h_trust; h_trust_x_slack1; h_trust_x_slack2; ...
        h_vc_slack1; h_vc_slack2; h_tau_slack1; h_tau_slack2; ...
        h_upperb_bounds; h_lowerb_bounds; h_lowerb_x0; h_upperb_x0; ...
        h_lowerb_xf; h_upperb_xf; h_lowerb_tau; ...
        h_upperb_tau; h_slack_ineq_bounds; h_socc_control];
    
    A_ecos = sparse(A_dynamics);
    b_ecos = b_dynamics;
    
    % Trust region (inequality constraints)
    % || xk - xk_ref ||_1 <= rk
    % this is the part of the trust region that changes -> we need to add the reference x to the h vector
    % and adjust the trust region radius
    for i=1:N % we consider all points
        xk = x_old(i, 1:n);
        h_trust_x_slack1(n*(i-1)+1: n*i) = -xk;
        h_trust_x_slack2(n*(i-1)+1: n*i) = xk;
        
        h_trust(i) = r_tr;
    end
    
    % Now we allocate those values in the final h vector
    h_ecos(h_trust_indices) = h_trust;
    h_ecos(h_trust_x_slack1_indices) = h_trust_x_slack1;
    h_ecos(h_trust_x_slack2_indices) = h_trust_x_slack2;
    
    % Call ECOS
    
    [ecos_result,~,info,~,~] = ecos(objective_c, G_ecos, h_ecos, dimensions_cones, A_ecos, b_ecos, opts);
    sum(abs(A_ecos*ecos_result - b_ecos))
    
    time = [time info.timing.runtime];
    
    % Post-processing
    virtual_ctrl = (reshape(ecos_result(x_len+u_len+1:x_len+u_len+virtual_ctrl_len),n,Ni*nc + 2));
    virtual_tau = ecos_result(x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+1: x_len+u_len+virtual_ctrl_len+aux_virtual_ctrl_len+virtual_tau_len);
    
    % Performance index
    J0_vect = zeros(1, sol_len);
    
    obj_coll_points = sum(auxdata.trans.weights_c_matrix*auxdata.phi.PHIu_matrix*auxdata.trans.Tu,1);
    obj_coll_points = obj_coll_points(x_len+1:x_len+u_len)';
    obj_node_points = sum(auxdata.trans.weights_n_matrix*auxdata.trans.Tu,1);
    obj_node_points = obj_node_points(x_len+1:x_len+u_len)';
    
    J0_vect(x_len+1:x_len+u_len) = h*0.5*(obj_coll_points+obj_node_points);
    
    J0 = J0_vect*ecos_result

    [paraECOS, paraGL, paraSCP, auxdata, hc_v, gc_v] = get_nonlinear_constr_viol(ecos_result, x_old, paraECOS, paraGL, paraSCP, auxdata);
    hc_v(1 : virtual_ctrl_len);
    
    % Sequential Convex Programming algorithm
    % Predicted objective function change
    phi_hat = J0 + lambda*sum(max(0, hc_v(virtual_ctrl_len + 1 : ...
        virtual_ctrl_len + virtual_tau_len))) + ratio_muc*muc*norm(virtual_ctrl,1) + ...
        ratio_lambda*lambda*sum(max(0,virtual_tau));
    
    % Actual objective function change
    phi_obj = J0 + ratio_muc*muc*norm(hc_v,1) + ratio_lambda*lambda*sum(max(0,gc_v));
    
    dphi = phi_obj_old - phi_obj
    dphi_hat = phi_hat_old - phi_hat;
    
     % Maximum constraint violation
    paraSCP.cmax = norm([hc_v, max(0,gc_v)]',Inf);
    cmax = paraSCP.cmax
    
    % Changing muc and lambda when feasibility is reached
    if cmax < epsc
        muc_n = muc/gamma
        lambda_n = lambda/gamma;
        ratio_muc = muc_n/muc;
        ratio_lambda = lambda_n/lambda;
    end

    % Rho definition
    rho = dphi/dphi_hat;
    
    if dphi < 0 && dphi_hat < 0
        rho = -1;
        dphi = 1;
    elseif dphi > 0 && dphi_hat < 0
        rho = rho1;
    elseif dphi < 0 && dphi_hat > 0
        rho = -1;
        dphi = 1;
    end
    
    % Update of alpha and beta: adaptive trust region radius update
    if k > 0
        if rho >= rho0 && rho_old >= rho0 % Case 1
            if beta < beta_max
                beta = beta*delta;
            end
            if alpha > delta
                alpha = alpha/delta;
            end
        elseif rho >= rho0 && rho_old < rho0 % Case 2
            if beta > delta
                beta = beta/delta;
            end
            if alpha < alpha_max
                alpha = alpha*delta;
            end
        elseif rho < rho0 && rho_old >= rho0 % Case 3
        elseif rho < rho0 && rho_old < rho0 % Case 4: Both rejected
            if alpha < alpha_max
                alpha = alpha*delta;
            end
        end
    end
    
    % Solution rejection/update
    if rho < rho0
        
        % Trust region radius
        r_tr = r_tr/alpha;
          
    else
        
        % Actual and predicted changes update
        
        phi_obj_old = phi_obj;
        phi_hat_old = phi_hat;
        
        % State update
        x_old_out = x_old;
        
        x_state_old = (reshape(ecos_result(1:x_len), n, N))';
        u_old = (reshape(ecos_result(x_len+1:x_len+u_len), m, N))';
        x_old = [x_state_old u_old];
        x = x_old;
        
        % Trust region radius update
        if rho < rho1
            r_tr = r_tr/alpha;
        elseif rho >= rho1 && rho < rho2
        elseif rho >= rho2
            r_tr = beta*r_tr;
        end
        
    end
    
    % Iterations
    k = k + 1
    
    % Alpha and beta vectors
    paraSCP.alpha_vect = [paraSCP.alpha_vect alpha];
    paraSCP.beta_vect = [paraSCP.beta_vect beta];
    
    % Update of rho
    rho_old = rho;
    
    % If feasibility is reached in iter_max iterations, the solution is
    % accepted
    if k == paraSCP.iter_max && cmax > epsc
        fprintf('Not converged');
        break
    elseif k == paraSCP.iter_max && cmax <= epsc
       dphi = epsphi*0.5;
    end
    
end
    
end


