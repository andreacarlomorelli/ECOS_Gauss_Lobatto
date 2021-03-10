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
th_f = auxdata.th_f;
th_0 = auxdata.th_0;
ni = auxdata.ni;
nf = auxdata.nf;
nrev = auxdata.nrev;

% GL parameters
np = paraGL.np;
Ni = paraGL.Ni;
N = paraGL.N;

Nseg = paraSCP.Nseg;

paraTRAJ.th_vect(:,1) = linspace(th_0,th_f + 2*nrev*pi,N); % Dimensionless time vector
paraTRAJ.th_vect_aux(:,1) = linspace(th_0,th_f + 2*nrev*pi, Ni + 1); % Auxiliary time vector for time step definition
h = paraTRAJ.th_vect_aux(2,1) - paraTRAJ.th_vect_aux(1,1); % Time step

% Initialization of the output parameters
paraTRAJ.Ni_iter = zeros(Nseg,1);
paraTRAJ.nr_iter = zeros(Nseg,1);

% SCP initialisation parameters
paraSCP.x0 = zeros(Nseg, ni); paraSCP.x0(1, :) = auxdata.x0;
paraSCP.xf = zeros(Nseg, nf); paraSCP.xf(1, :) = auxdata.xf;
paraTRAJ.nr_iter(1) = round(N/Nseg); paraTRAJ.sw_nodes = zeros(Nseg,1); 
th_iter = th_0; paraTRAJ.h_vect = h; h_iter = h; th_f_iter = th_f + 2*nrev*pi; 
paraSCP.th_f_final = th_f;

% Initial guess phi and phi_hat
[paraECOS, paraGL, paraSCP, auxdata] = get_phi_initial_guess(paraECOS, paraGL, paraSCP, auxdata, paraSCP.x_old);

for e = 1 : Nseg
    
    paraSCP.e = e;
    
    % Initialisation parameters
    paraSCP.r_tr = paraSCP.r0*e;
    dth = (th_f + 2*nrev*pi)/(100*e);

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
        
        th_iter = th_iter + h_iter*paraTRAJ.Ni_iter(e);
        
        paraSCP.x_old = interp1(paraTRAJ.th_vect(:,e), x(:,:), ...
            linspace(th_iter, th_f_iter, N));

        % Propagation of the planets' state
        [paraTRAJ.x_cart(:,:,e)] = transform_state_cyl_to_cart(paraTRAJ.x(:,:,e), linspace(th_iter, th_f_iter,N));
        [~,state_planet] = ode113(@f_planet, [th_f_iter th_f_iter + dth], paraTRAJ.x_cart(end,:,e));
        [state_planet_cyl] = transform_state_cart_to_cyl(state_planet);
        
        th_f_iter = th_f_iter + dth;

        paraSCP.th_f_final = th_f_iter;

        paraSCP.xf(e+1,1) = state_planet_cyl(end,1); paraSCP.xf(e+1,2) = state_planet_cyl(end,2);
        paraSCP.xf(e+1,3) = state_planet_cyl(end,3); paraSCP.xf(e+1,4) = state_planet_cyl(end,4);
        paraSCP.xf(e+1,5) = state_planet_cyl(end,5);
        
        % (Perturbed) initial BCs for next step
        paraSCP.x0(e + 1, 1 : 2) = paraSCP.x_old(1, 1 : 2) + (1e8/auxdata.units.R0) - (2e8/auxdata.units.R0)*rand(1, 2);
        paraSCP.x0(e + 1, 3 : 5) = paraSCP.x_old(1, 3 : 5) + (1e3/auxdata.units.V0) - (2e3/auxdata.units.R0)*rand(1, 3);
        paraSCP.x0(e + 1, 6) = paraSCP.x_old(1, 6);
        paraSCP.x0(e + 1, 7) = paraSCP.x_old(1, 7);

        paraTRAJ.th_vect(:,e + 1) = linspace(th_iter, th_f_iter, N);
        paraTRAJ.th_vect_aux(:, e + 1) = linspace(th_iter, th_f_iter, Ni + 1);
        h_iter = paraTRAJ.th_vect_aux(2, e + 1) - ...
            paraTRAJ.th_vect_aux(1, e + 1);
        paraTRAJ.h_vect(e + 1) = h_iter;
        paraSCP.h = h_iter;
              
    end
    
end

end

