function [paraGL, paraTRAJ, paraSCP, auxdata] = dynrec(paraGL, paraTRAJ, paraSCP, auxdata)
%
% dynrec
% 
% FUNCTION DESCRIPTION: this function reconstructs the dynamics of the
% S/C with the given thrust hustory.
% 
% INPUTS:
% 
% paraTRAJ:    structure containing parameters and variables related to
%              the trajectory
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

% Sates and controls
n = 7;
m = 4;

% GL parameters
nc = paraGL.nc;
Ni = paraGL.Ni;

% auxdata parameters
c = auxdata.engine.c;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

% SCP parameters
Nseg = paraSCP.Nseg;

options = odeset('AbsTol',1e-12, 'RelTol', 1e-12);
paraTRAJ.x_rec = zeros(nc*Ni + 2, n + m);

Ni_tot = sum(paraTRAJ.Ni_iter);
paraTRAJ.x_g_rec = zeros(nc*Ni_tot + 2 + (Nseg - 1), n);
paraTRAJ.x_rec = zeros(nc*Ni + 2, n, Nseg);

for e = 1 : Nseg % Loop for all the trajectory segments
    
    % Initial guess for trajectory reconstruction
    x0_rec = [paraTRAJ.x_c(1,1,e) paraTRAJ.x_c(1,2,e) ...
        paraTRAJ.x_c(1,3,e) paraTRAJ.x_c(1,4,e) paraTRAJ.x_c(1,5,e) ...
        paraTRAJ.x_c(1,6,e) paraTRAJ.x_c(1,7,e)];
    
    taux = paraTRAJ.u_c(:,1,e); tauy = paraTRAJ.u_c(:,2,e);
    tauw = paraTRAJ.u_c(:,3,e); tau = paraTRAJ.u_c(:,4,e);
    
    % ode113
    [~,paraTRAJ.x_rec(:,:,e)] = ode113(@dyn_rec, paraTRAJ.t_vect_coll(:,e), x0_rec, ...
        options, c, ve, V0, paraTRAJ.t_vect_coll(:,e), tau, taux, tauy, tauw);
    
    % Global trajectory assembly
    if Nseg == 1
        paraTRAJ.x_g_rec(1 : nc*paraTRAJ.Ni_iter(e) + 2, :) = ...
            paraTRAJ.x_rec(1 : nc*paraTRAJ.Ni_iter(e) + 2,:,e);
    elseif e == Nseg
        paraTRAJ.x_g_rec(a + 1: a + 1 + nc*paraTRAJ.Ni_iter(e) + 1, :) = ...
            paraTRAJ.x_rec(1 : nc*paraTRAJ.Ni_iter(e) + 2,:,e);
        a = a + 1 + nc*paraTRAJ.Ni_iter(e) + 1;
    elseif e == 1
        paraTRAJ.x_g_rec(1 : nc*paraTRAJ.Ni_iter(1) + 1, :) = ...
            paraTRAJ.x_rec(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e);
        a =  nc*paraTRAJ.Ni_iter(1) + 1;
    else
        paraTRAJ.x_g_rec(a + 1 : a + 1 + nc*paraTRAJ.Ni_iter(e), :) = ...
            paraTRAJ.x_rec(1 : nc*paraTRAJ.Ni_iter(e) + 1,:,e);
        a = a + 1 + nc*paraTRAJ.Ni_iter(e);
    end
    
end

end

