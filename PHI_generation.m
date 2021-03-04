function [paraGL, auxdata] = PHI_generation(paraGL, auxdata)
%
% PHI_generation
% 
% FUNCTION DESCRIPTION: this function computes the PHI matrices associated
% to the Gauss-Lobatto method.
% 
% INPUTS/OUTPUTS: 
%
% paraGL:        Structure with the GL parameters  
%
% auxdata:       Structure with the auxiliary parameters
% 
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 1/03/2021
%

% Sates and controls
n = 7;
m = 4;

% GL parameters
nc = paraGL.nc;
np = paraGL.np;
Ni = paraGL.Ni;

% Baseline matrices
PHI_c = paraGL.PHI_c;
PHI_pc = paraGL.PHI_pc;
PHI_n = paraGL.PHI_n; 
PHI_pn = paraGL.PHI_pn; 
PHI_u = paraGL.PHI_u;
weights_c = paraGL.weights_c';
weights_n = paraGL.weights_n';

% Initializing the matrices
auxdata.phi.PHI_matrix = zeros(nc*n*Ni, 2*np*Ni*n);
auxdata.phi.PHI_p_matrix = zeros(nc*n*Ni, 2*np*Ni*n);
auxdata.phi.PHIu_matrix = zeros(nc*m*Ni, np*Ni*m);
auxdata.phi.PHIpn_matrix_1 = zeros(n, np*n);
auxdata.phi.PHIpn_matrix_end = zeros(n, np*n);
auxdata.phi.PHIn_matrix_1 = zeros(n, np*n);
auxdata.phi.PHIn_matrix_end = zeros(n, np*n);
auxdata.phi.PHIun_matrix_1 = zeros(m, np*m);
auxdata.trans.weights_c_matrix = zeros(Ni, nc*Ni*m);
auxdata.trans.weights_n_matrix = zeros(Ni, nc*Ni*m);

% PHI_matrix
PHI_big = zeros(n*nc, 2*np*n);
for j = 1 : nc
    for i = 1 : n
        PHI_big((j-1)*n + i, (i-1)*2*np + 1 : 2*np*i) = PHI_c(j,:);
    end
end
PHI_cell = repmat({PHI_big},1,Ni);
PHI_matrix_par = blkdiag(PHI_cell{:});
auxdata.phi.PHI_matrix(:, 1 : size(PHI_matrix_par,2)) = PHI_matrix_par;

% PHI_p_matrix
PHI_p_big = zeros(n*nc, 2*np*n);
for j = 1 : nc
    for i = 1 : n
        PHI_p_big((j-1)*n + i, (i-1)*2*np + 1 : 2*np*i) = PHI_pc(j,:);
    end
end
PHI_p_cell = repmat({PHI_p_big},1,Ni);
PHI_p_matrix_par = blkdiag(PHI_p_cell{:});
auxdata.phi.PHI_p_matrix(:, 1 : size(PHI_p_matrix_par,2)) = PHI_p_matrix_par;

% PHIu_matrix
PHIu_big = zeros(m*nc, np*m);
for j = 1 : nc
    for i = 1 : m
        PHIu_big((j-1)*m + i, (i-1)*np + 1 : np*i) = PHI_u(j,:);
    end
end
PHIu_cell = repmat({PHIu_big},1, Ni);
PHIu_matrix_par = blkdiag(PHIu_cell{:});
auxdata.phi.PHIu_matrix(:, 1 : size(PHIu_matrix_par,2)) = PHIu_matrix_par;

% PHIn_matrix_1
PHIpn_cell_1 = repmat({PHI_pn(1,:)},1,n);
PHIpn_matrix_par_1 = blkdiag(PHIpn_cell_1{:});
auxdata.phi.PHIpn_matrix_1(:, 1 : size(PHIpn_matrix_par_1,2)) = PHIpn_matrix_par_1;

PHIn_cell_1 = repmat({PHI_n(1,:)},1,n);
PHIn_matrix_par_1 = blkdiag(PHIn_cell_1{:});
auxdata.phi.PHIn_matrix_1(:, 1 : size(PHIn_matrix_par_1,2)) = PHIn_matrix_par_1;

% PHIn_matrix_end
PHIpn_cell_end = repmat({PHI_pn(end,:)},1,n);
PHIpn_matrix_par_end = blkdiag(PHIpn_cell_end{:});
auxdata.phi.PHIpn_matrix_end(:, 1 : size(PHIpn_matrix_par_end,2)) = PHIpn_matrix_par_end;

PHIn_cell_end = repmat({PHI_n(end,:)},1,n);
PHIn_matrix_par_end = blkdiag(PHIn_cell_end{:});
auxdata.phi.PHIn_matrix_end(:, 1 : size(PHIn_matrix_par_end,2)) = PHIn_matrix_par_end;

% PHIun_matrix_1
PHIun_cell_1 = repmat({PHI_u(1,:)},1,m);
PHIun_matrix_par_1 = blkdiag(PHIun_cell_1{:});
auxdata.phi.PHIun_matrix_1(:, 1 : size(PHIun_matrix_par_1,2)) = PHIun_matrix_par_1;

% PHIun_matrix_end
PHIun_cell_end = repmat({PHI_u(end,:)},1,m);
PHIun_matrix_par_end = blkdiag(PHIun_cell_end{:});
auxdata.phi.PHIun_matrix_end(:, 1 : size(PHIun_matrix_par_end,2)) = PHIun_matrix_par_end;

% Weights for objective function
weights_c_cell_par = zeros(1,m*nc);
for i = 1 : nc
    weights_c_cell_par(m*i) = weights_c(i);
end

weights_c_cell = repmat({weights_c_cell_par},1,Ni);
auxdata.trans.weights_c_matrix =  blkdiag(weights_c_cell{:});

weights_n_cell_par = zeros(1,np*m);
weights_n_cell_par((m-1)*np + 1 : end) = weights_n;

weights_n_cell = repmat({weights_n_cell_par},1,Ni);
auxdata.trans.weights_n_matrix =  blkdiag(weights_n_cell{:});

end

