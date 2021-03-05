%% Convex Optimization fuel-optimal transfer with Gauss - Lobatto quadrature rule and ECOS
clc; clear all; 

% This script implements the Legendre-Gauss-Lobatto discretization scheme
% applied to the low-thrust convex optimization problem and uses the ECOS
% solver to solve it.

% Test

%% Problem parameters

% Transfer-related data load from function
transfer = 1;
if transfer == 1
    auxdata = Earth_Mars_parameters;
elseif transfer == 2
    auxdata = Earth_Dionysus_parameters;
elseif transfer == 3
    auxdata = Earth_Venus_parameters;
end

% Variables bounds-related data load from function
auxdata = bounds(auxdata);

% Number of trajectory segments
paraSCP.Nseg = 5;

% Gauss-Lobatto method order
paraGL.ng = input('Number of Gauss - Lobatto points: ');

% Number of segments
paraGL.Ni = input('Number of segments ');

% Initial and final BCs
t0 = auxdata.t0;
tf = auxdata.tf;
rr_0 = auxdata.rr_0;
vv_0 = auxdata.vv_0;
rr_f = auxdata.rr_f;
vv_f = auxdata.vv_f;

% Initial and final conditions vectors
auxdata.x0 = [rr_0(1) rr_0(2) rr_0(3) vv_0(1) vv_0(2) vv_0(3) 0];
auxdata.xf = [rr_f(1) rr_f(2) rr_f(3) vv_f(1) vv_f(2) vv_f(3)];
paraSCP.x0 = auxdata.x0;

auxdata.ni = length(auxdata.x0);
auxdata.nf = length(auxdata.xf);

% number of states (n = 7) and controls (m = 4)
n = 7;
m = 4;

% Useful parameters
c = auxdata.engine.c;
ve = auxdata.engine.ve;

% Gauss - Lobatto parameters & matrices
paraGL = GL(paraGL);

np = paraGL.np;
nc = paraGL.nc;
N = paraGL.N;
Ni = paraGL.Ni;

% Generation of the PHI matrices
[paraGL, auxdata] = PHI_generation(paraGL, auxdata);

% Time vectors
t_vect_adim = linspace(0,tf,N); % Dimensionless time vector
t_vect_aux = linspace(0,tf, Ni + 1); % Auxiliary time vector for time step definition

% Time step
paraSCP.h = t_vect_aux(2) - t_vect_aux(1);
h = paraSCP.h;

%% Initialization of the SCP algorithm parameters

% Penalty weights
if transfer == 2
    paraSCP.lambda = 500;
    paraSCP.muc = 500;
else
    paraSCP.lambda = 100;
    paraSCP.muc = 100;
end

% Initial trust region radius
paraSCP.r0 = input('Initial trust region radius ');

paraSCP.rho0 = 0.01;
paraSCP.rho1 = 0.25;
paraSCP.rho2 = 0.9;

paraSCP.alpha = 1.4;
paraSCP.alpha_max = 4;
paraSCP.beta = 1.4;
paraSCP.beta_max = 4;
paraSCP.delta = 1.3;
paraSCP.gamma = 1;

% Convergence thresholds
if transfer == 2
    paraSCP.epsc = 1e-6;
    paraSCP.epsphi = 1e-2;
else
    paraSCP.epsc = 1e-6;
    paraSCP.epsphi = 1e-4;
end

% Initial guess
inguess = 2; % 1 for Cubic - based, 2 for FFS
nrev = input('How many revolutions do you want to consider for your initial guess? ');
auxdata.Ta_max = 100; % Maximum thrust for FFS initial guess

[x_guess] = inguess_ecos(N, nrev, inguess, auxdata);
paraSCP.x_old = x_guess;

% Maximum number of iterations
paraSCP.iter_max = 500;

%% ECOS transformation matrices

% Solution vector's length definition
% The vector is:
% [states, controls, virtual controls, ...
% auxiliary variables for virtual controls, ...
% virtual tau, auxiliary variables for virtual tau, ...
% auxiliary variables for trust region]

% State variables
x_len = n*N;

% Control variables
u_len = m*N;

% State variables at collocation points (not including first and last)
x_len_c = nc*Ni*n;

% Control variables at collocation points (not including first and last)
u_len_c = nc*Ni*m;

% Virtual variables for the inequality constraint on tau
virtual_tau_len = nc*Ni + 2;
aux_virtual_tau_len = virtual_tau_len;

% Virtual control variables 
virtual_ctrl_len = nc*Ni*n +2*n;
aux_virtual_ctrl_len = virtual_ctrl_len;

% Total number of virtual variables
slack_ineq_len = virtual_ctrl_len + aux_virtual_ctrl_len + virtual_tau_len + aux_virtual_tau_len;

% Trust region radius variables
aux_trust_x_len = x_len;

% Total length of the vector: 
sol_len = x_len + u_len + virtual_ctrl_len + aux_virtual_ctrl_len + virtual_tau_len + aux_virtual_tau_len + aux_trust_x_len;

% Vector of lengths
paraECOS.len_vect = [x_len, u_len, x_len_c, u_len_c, virtual_tau_len, virtual_ctrl_len, aux_virtual_ctrl_len, aux_virtual_tau_len, slack_ineq_len, aux_trust_x_len, sol_len];

% Generation of the constant parts of the martices T and Tu
[paraECOS, paraGL, auxdata] = get_constant_T_Tu(paraECOS, paraGL, auxdata);

% Generation of the varying parts of the matrices T and Tu
[paraECOS, paraGL, paraSCP, auxdata] = get_varying_T_Tu(x_guess, paraECOS, paraGL, paraSCP, auxdata);

% Generation of the constant parts of the matrix G
paraSCP.e = 1;
[paraECOS, paraGL, paraSCP, auxdata] = get_constant_matrices(paraECOS, paraGL, paraSCP, auxdata);

% Get varying parts of matrices A_ecos, G_ecos, h_ecos and b_ecos
[paraGL, paraSCP, paraECOS, auxdata] = get_varying_matrices(x_guess, paraGL, paraSCP, paraECOS, auxdata);

%% ECOS Sequential Convex Programming

tic
[paraECOS, paraSCP, paraTRAJ, paraGL, auxdata] = trajectory(paraECOS, paraSCP, paraGL, auxdata);
toc

%% Post-processing

% Assembling the different trajectory segments at the nodes and collocation
% points
[paraGL, paraSCP, paraTRAJ, auxdata] = postproc(paraGL, paraSCP, paraTRAJ, auxdata);

% Time vector at the collocation points
t_g = paraTRAJ.t_g;

% Time vector in days
R0 = auxdata.units.R0;
V0 = auxdata.units.V0;

t_g_days = t_g*(R0/V0)/(24*3600);

% Calling variables at collocation points
xvar_c = paraTRAJ.x_g(:,1); y_c = paraTRAJ.x_g(:,2); w_c = paraTRAJ.x_g(:,3);
vx_c = paraTRAJ.x_g(:,4); vy_c = paraTRAJ.x_g(:,5); vw_c = paraTRAJ.x_g(:,6);
z_c = paraTRAJ.x_g(:,7); taux_c = paraTRAJ.x_g(:,8); tauy_c = paraTRAJ.x_g(:,9);
tauw_c = paraTRAJ.x_g(:,10); tau_c = paraTRAJ.x_g(:,11);

% Mass
mass = exp(z_c)*auxdata.sc.m0;

figure
plot3(xvar_c, y_c, w_c)

% Dynamics reconstruction
[paraGL, paraTRAJ, paraSCP, auxdata] = dynrec(paraGL, paraTRAJ, paraSCP, auxdata);

% Calling reconstructed variables at collocation points
xvar_rec = paraTRAJ.x_g_rec(:,1); y_rec = paraTRAJ.x_g_rec(:,2); 
w_rec = paraTRAJ.x_g_rec(:,3); vx_rec = paraTRAJ.x_g_rec(:,4); 
vy_rec = paraTRAJ.x_g_rec(:,5); vw_rec = paraTRAJ.x_g_rec(:,6); 
z_rec = paraTRAJ.x_g_rec(:,7);

% Variabvles errors as functions of ToF
err_x = abs(xvar_c - xvar_rec);
err_y = abs(y_c - y_rec);
err_w = abs(w_c - w_rec);
err_vx = abs(vx_c - vx_rec);
err_vy = abs(vy_c - vy_rec);
err_vw = abs(vw_c - vw_rec);
err_z = abs(z_c - z_rec);

% Errors w.r.t the final boundary conditions
x_bound = abs(xvar_rec(end) - auxdata.xf(1));
y_bound = abs(y_rec(end) - auxdata.xf(2));
w_bound = abs(w_rec(end) - auxdata.xf(3));
vx_bound = abs(vx_rec(end) - auxdata.xf(4));
vy_bound = abs(vy_rec(end) - auxdata.xf(5));
vw_bound = abs(vw_rec(end) - auxdata.xf(6));

% Average boundary errors
bound_pos = mean([x_bound, y_bound, w_bound])
bound_vel = mean([vx_bound, vy_bound, vw_bound])

%% Thrust plot

% Thrust of the complete trajectory
figure
T_c = exp(paraTRAJ.x_g(:,7)).*paraTRAJ.x_g(:,11);
hold all
plot(t_g_days,T_c,'k','LineWidth',2)
grid on

figure
hold all
plot(taux_c)
plot(tauy_c)
plot(tauw_c)
plot(tau_c)
 