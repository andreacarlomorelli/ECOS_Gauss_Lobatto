function [paraGL, paraSCP, paraTRAJ, auxdata] = postproc(paraGL, paraSCP, paraTRAJ, auxdata)
%
% postproc
% 
% FUNCTION DESCRIPTION: this function assembles the various trajectory 
% segments.
%  
% INPUTS:
% 
% paraTRAJ:    structure containing parameters and variables related to
%              the trajectory
%
% paraGL:      structure containing parameters and variables related to
%              the Gauss-Lobatto method.
%
% Nseg:        Number of trajectory segments [1x1]
%
% OUTPUTS:
%
% paraTRAJ:    structure containing parameters and variables related to
%              the trajectory
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 26/11/2020
%

% number of states (n = 7) and controls (m = 4)
n = 7;
m = 4;

% auxdata parameters
th_f = auxdata.th_f;
nrev = auxdata.nrev;

% GL parameters
nc = paraGL.nc;
np = paraGL.np;
Ni = paraGL.Ni;
N = paraGL.N;

% SCP parameters
Nseg = paraSCP.Nseg;

% Time
paraTRAJ.t_vect_coll = zeros(Ni*nc + 2, Nseg);

% Variables at nodes
paraTRAJ.xvar = zeros(N, Nseg); paraTRAJ.y = zeros(N, Nseg);
paraTRAJ.w = zeros(N, Nseg); paraTRAJ.vx = zeros(N, Nseg);
paraTRAJ.vy = zeros(N, Nseg); paraTRAJ.vw = zeros(N, Nseg);
paraTRAJ.z = zeros(N, Nseg); paraTRAJ.taux = zeros(N, Nseg);
paraTRAJ.tauy = zeros(N, Nseg); paraTRAJ.tauw = zeros(N, Nseg);
paraTRAJ.tau = zeros(N, Nseg);

% Collocation variables
paraTRAJ.xvar_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.y_c = zeros(nc*Ni + 2, Nseg);
paraTRAJ.w_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.vx_c = zeros(nc*Ni + 2, Nseg);
paraTRAJ.vy_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.vw_c = zeros(nc*Ni + 2, Nseg);
paraTRAJ.z_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.taux_c = zeros(nc*Ni + 2, Nseg);
paraTRAJ.tauy_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.tauw_c = zeros(nc*Ni + 2, Nseg);
paraTRAJ.tau_c = zeros(nc*Ni + 2, Nseg); paraTRAJ.T = zeros(nc*Ni + 2, Nseg);

paraTRAJ.x_c = zeros(nc*Ni + 2, n,Nseg);
paraTRAJ.u_c = zeros(nc*Ni + 2, m,Nseg);

% Global trajectory
Ni_tot = sum(paraTRAJ.Ni_iter);
paraTRAJ.x_g = zeros(nc*Ni_tot + 2 + (Nseg - 1), n + m);
paraTRAJ.x_g(end, :) = paraTRAJ.x(end, :,end);
paraTRAJ.t_g = zeros(nc*Ni_tot + 2 + (Nseg - 1),1);
paraTRAJ.t_g(1,1) = 0;
paraTRAJ.t_g(end,1) = th_f + 2*nrev*pi;

for e = 1 : Nseg
    
    h_rec = paraTRAJ.h_vect(e);
    
    % Time vector
    paraTRAJ.t_vect_coll(1,e) = paraTRAJ.t_vect(1,e);
    paraTRAJ.t_vect_coll(end,e) =  paraTRAJ.t_vect(end,e);
    
    for j = 1 : Ni
        paraTRAJ.t_vect_coll(2 + (j-1)*nc : 1 + j*nc, e) = paraGL.zeta*0.5*h_rec + ...
            (paraTRAJ.t_vect_aux(j + 1,e) + paraTRAJ.t_vect_aux(j,e))*0.5;
    end
    
    % Nodal variables
    paraTRAJ.xvar(:,e) = paraTRAJ.x(:,1,e); paraTRAJ.y(:,e) = paraTRAJ.x(:,2,e);
    paraTRAJ.w(:,e) = paraTRAJ.x(:,3,e); paraTRAJ.vx(:,e) = paraTRAJ.x(:,4,e);
    paraTRAJ.vy(:,e) = paraTRAJ.x(:,5,e); paraTRAJ.vw(:,e) = paraTRAJ.x(:,6,e);
    paraTRAJ.z(:,e) = paraTRAJ.x(:,7,e); paraTRAJ.taux(:,e) = paraTRAJ.x(:,8,e);
    paraTRAJ.tauy(:,e) = paraTRAJ.x(:,9,e); paraTRAJ.tauw(:,e) = paraTRAJ.x(:,10,e);
    paraTRAJ.tau(:,e)= paraTRAJ.x(:,11,e); paraTRAJ.T = zeros(nc*Ni + 2, Nseg);
    
    for j = 1 : Ni
        
        [b_c_x, b_tau] = collpoints(np, h_rec, paraTRAJ.x(:,:,e), paraTRAJ.x_old(:,:,e), j, auxdata);
        
        % State
        paraTRAJ.x_c((j-1)*nc + 2 : j*nc + 1,:,e) = paraGL.PHI_c*b_c_x;
        
        % Control
        paraTRAJ.u_c((j-1)*nc + 2 : j*nc + 1,:,e) = paraGL.PHI_u*b_tau;
        
    end
    
    paraTRAJ.x_c(1,:,e) = paraTRAJ.x(1, 1 : n, e);
    paraTRAJ.x_c(end,:,e) = paraTRAJ.x(end, 1 : n, e);
    
    paraTRAJ.u_c(1,:,e) = paraTRAJ.x(1, n + 1 : n + m, e);
    paraTRAJ.u_c(end,:,e) = paraTRAJ.x(end, n + 1 : n + m, e);
    
    paraTRAJ.xvar_c(:,e) = paraTRAJ.x_c(:,1,e); paraTRAJ.y_c(:,e) = paraTRAJ.x_c(:,2,e);
    paraTRAJ.w_c(:,e) = paraTRAJ.x_c(:,3,e); paraTRAJ.vx_c(:,e) = paraTRAJ.x_c(:,4,e);
    paraTRAJ.vy_c(:,e) = paraTRAJ.x_c(:,5,e); paraTRAJ.vw_c(:,e) = paraTRAJ.x_c(:,6,e);
    paraTRAJ.z_c(:,e) = paraTRAJ.x_c(:,7,e);
    
    paraTRAJ.taux_c(:,e) = paraTRAJ.u_c(:,1,e); paraTRAJ.tauy_c(:,e) = paraTRAJ.u_c(:,2,e);
    paraTRAJ.tauw_c(:,e) = paraTRAJ.u_c(:,3,e); paraTRAJ.tau_c(:,e) = paraTRAJ.u_c(:,4,e);
    
    paraTRAJ.T(:,e) = paraTRAJ.tau_c(:,e).*exp(paraTRAJ.z_c(:,e)); % Non dimensional thrust
    
    % Global trajectory assembly
    if e == 1
        paraTRAJ.x_g(1 : nc*paraTRAJ.Ni_iter(e) + 1, :) = ...
            [paraTRAJ.x_c(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e) ...
            paraTRAJ.u_c(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e)];
        paraTRAJ.t_g(1 : nc*paraTRAJ.Ni_iter(e) + 1, :) = ...
            paraTRAJ.t_vect_coll(1 : nc*paraTRAJ.Ni_iter(e) + 1,e);
       
        a = nc*paraTRAJ.Ni_iter(e) + 1;
    else
        paraTRAJ.x_g(a + 1 : ...
            a + 1 + nc*paraTRAJ.Ni_iter(e), :) = ...
            [paraTRAJ.x_c(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e) ...
            paraTRAJ.u_c(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e)];
        paraTRAJ.t_g(a + 1 : a + 1 + ...
            nc*paraTRAJ.Ni_iter(e)) = paraTRAJ.t_vect_coll...
            (1 : nc*paraTRAJ.Ni_iter(e) + 1,e);
        
        a = a + 1 + nc*paraTRAJ.Ni_iter(e);
    end
end

end

