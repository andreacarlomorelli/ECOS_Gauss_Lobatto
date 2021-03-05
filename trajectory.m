function [paraECOS, paraSCP, paraTRAJ, paraGL, auxdata] = trajectory(paraECOS, paraSCP, paraGL, auxdata)
%
% trajectory
% 
% FUNCTION DESCRIPTION: this function applies the SCP method to
% each of the trajectory segments.
% 
% INPUTS:
% 
% paraSCP      Structure containing parameters and variables related to the
%              Sequential Convex Programming algorithm
%
% paraGL       Structure containing parameters and variables related to the
%              Gauss-Lobatto method.
%
% x_guess      Vector of the initial trajectory guess   [m_guess x 11]
%
% OUTPUTS:
%
% paraSCP      Structure containing parameters and variables related to the
%              Sequential Convex Programming algorithm
%
% paraGL       Structure containing parameters and variables related to the
%              Gauss-Lobatto method.
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 27/11/2020
%

% auxdata parameters
tf = auxdata.tf;
ni = auxdata.ni;
nf = auxdata.nf;

% GL parameters
np = paraGL.np;
Ni = paraGL.Ni;
N = paraGL.N;

Nseg = paraSCP.Nseg;

paraTRAJ.t_vect(:,1) = linspace(0,tf,N); % Dimensionless time vector
paraTRAJ.t_vect_aux(:,1) = linspace(0,tf, Ni + 1); % Auxiliary time vector for time step definition
h = paraTRAJ.t_vect_aux(2,1) - paraTRAJ.t_vect_aux(1,1); % Time step

% Initialization of the output parameters
paraTRAJ.Ni_iter = zeros(Nseg,1);
paraTRAJ.nr_iter = zeros(Nseg,1);

% SCP initialisation parameters
paraSCP.x0 = zeros(Nseg, ni); paraSCP.x0(1, :) = auxdata.x0;
paraSCP.xf = zeros(Nseg, nf); paraSCP.xf(1, :) = auxdata.xf;
paraTRAJ.nr_iter(1) = round(N/Nseg); paraTRAJ.sw_nodes = zeros(Nseg,1); 
TOF_iter = 0; paraTRAJ.h_vect = h; h_iter = h; tf_iter = tf; 

% Initial guess phi and phi_hat
[paraECOS, paraGL, paraSCP, auxdata] = get_phi_initial_guess(paraECOS, paraGL, paraSCP, auxdata, paraSCP.x_old);

for e = 1 : Nseg
    
    paraSCP.e = e;
    
    % Initialisation parameters
    paraSCP.r_tr = paraSCP.r0*e;
    dtf = tf/(100*e);

    % Number of segments and nodes at iteration e
    paraTRAJ.Ni_iter(e) = round(Ni/(Nseg - e + 1));
    paraTRAJ.nr_iter(e) = np + (np - 1)*(paraTRAJ.Ni_iter(e) - 1);
    
    % SCP algorithm
    
    % Generation of the constant parts of the matrix G

    [paraECOS, paraGL, paraSCP, auxdata] = get_constant_matrices(paraECOS, paraGL, paraSCP, auxdata);
    
    [paraECOS, paraGL, paraSCP, x, x_old, time, J0] = ECOS_SCP(paraECOS, paraGL, paraSCP, auxdata);
    
    % Solution at iteration e
    paraTRAJ.x(:,:,e) = x;
    paraTRAJ.x_old(:,:,e) = x_old;

    % Switching nodes
    paraTRAJ.sw_nodes(e) = paraTRAJ.nr_iter(e);
    
    % Settings for next iteration
    if e < Nseg
        
        TOF_iter = TOF_iter + h_iter*paraTRAJ.Ni_iter(e);
        
        paraSCP.x_old = interp1(paraTRAJ.t_vect(:,e), x(:,:), ...
            linspace(TOF_iter, tf, N));

        % Propagation of the planets' state
        [~,state_planet] = ode113(@f_planet, [tf_iter tf_iter + dtf], paraSCP.xf(e,:)');

        tf_iter = tf_iter + dtf;

        paraSCP.tf_final = tf_iter;

        paraSCP.xf(e+1,1) = state_planet(end,1); paraSCP.xf(e+1,2) = state_planet(end,2);
        paraSCP.xf(e+1,3) = state_planet(end,3); paraSCP.xf(e+1,4) = state_planet(end,4);
        paraSCP.xf(e+1,5) = state_planet(end,5); paraSCP.xf(e+1,6) = state_planet(end,6);
     
        % (Perturbed) initial BCs for next step
        paraSCP.x0(e + 1, 1 : 3) = paraSCP.x_old(1, 1 : 3) + (1e8/auxdata.units.R0) - (2e8/auxdata.units.R0)*rand(1, 3);
        paraSCP.x0(e + 1, 4 : 6) = paraSCP.x_old(1, 4 : 6) + (1e3/auxdata.units.V0) - (2e3/auxdata.units.R0)*rand(1, 3);
        paraSCP.x0(e + 1, 7) = paraSCP.x_old(1, 7);

        paraTRAJ.t_vect(:,e + 1) = linspace(TOF_iter, tf_iter, N);
        paraTRAJ.t_vect_aux(:, e + 1) = linspace(TOF_iter, tf_iter, Ni + 1);
        h_iter = paraTRAJ.t_vect_aux(2, e + 1) - ...
            paraTRAJ.t_vect_aux(1, e + 1);
        paraTRAJ.h_vect(e + 1) = h_iter;
        paraSCP.h = h_iter;
              
    end
    
end

end

