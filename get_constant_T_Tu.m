function [paraECOS, paraGL, auxdata] = get_constant_T_Tu(paraECOS, paraGL, auxdata)
%
% get_constant_T_Tu
% 
% FUNCTION DESCRIPTION: this function computes the constant parts of the
% matrices T and Tu
% 
% INPUTS/OUTPUTS:
% 
% paraECOS:      Structure with the ECOS parameters  
%
% paraGL:        Structure with the GL parameters  
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
sol_len = len_vect(11);

% GL parameters
np = paraGL.np;
Ni = paraGL.Ni;
N = paraGL.N;

% T matrices
T = zeros(2*np*n*Ni, sol_len);
Tu = zeros(np*m*Ni, sol_len);

Taux_x = [1 0 0 0 0 0 0];
Tcell_x = repmat({Taux_x},1,np);
Taux_y = [0 1 0 0 0 0 0];
Tcell_y = repmat({Taux_y},1,np);
Taux_w = [0 0 1 0 0 0 0];
Tcell_w = repmat({Taux_w},1,np);
Taux_vx = [0 0 0 1 0 0 0];
Tcell_vx = repmat({Taux_vx},1,np);
Taux_vy = [0 0 0 0 1 0 0];
Tcell_vy = repmat({Taux_vy},1,np);
Taux_vw = [0 0 0 0 0 1 0];
Tcell_vw = repmat({Taux_vw},1,np);
Taux_z = [0 0 0 0 0 0 1];
Tcell_z = repmat({Taux_z},1,np);

Taux_taux = [1 0 0 0];
Tcell_taux = repmat({Taux_taux},1,np);
Taux_tauy = [0 1 0 0];
Tcell_tauy = repmat({Taux_tauy},1,np);
Taux_tauw = [0 0 1 0];
Tcell_tauw = repmat({Taux_tauw},1,np);
Taux_tau = [0 0 0 1];
Tcell_tau = repmat({Taux_tau},1,np);

for i = 1 : Ni
 
    % T matrix
    T((i-1)*2*np*n + 1 : (i-1)*2*np*n + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_x{:});
    T((i-1)*2*np*n + 2*np + 1 : (i-1)*2*np*n + 2*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_y{:});
    T((i-1)*2*np*n + 4*np + 1 : (i-1)*2*np*n + 4*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_w{:});
    T((i-1)*2*np*n + 6*np + 1 : (i-1)*2*np*n + 6*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_vx{:});
    T((i-1)*2*np*n + 8*np + 1 : (i-1)*2*np*n + 8*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_vy{:});
    T((i-1)*2*np*n + 10*np + 1 : (i-1)*2*np*n + 10*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_vw{:});
    T((i-1)*2*np*n + 12*np + 1 : (i-1)*2*np*n + 12*np + np, ...
        (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = blkdiag(Tcell_z{:});
    
    % Tu matrix
    Tu((i-1)*np*m + 1 : (i-1)*np*m + np, ...
        n*N + (i-1)*m*np - (i-1)*m + 1 : n*N + i*m*np - (i-1)*m) = blkdiag(Tcell_taux{:});
    Tu((i-1)*np*m + np + 1 : (i-1)*np*m + 2*np, ...
        n*N + (i-1)*m*np - (i-1)*m + 1 : n*N + i*m*np - (i-1)*m) = blkdiag(Tcell_tauy{:});
    Tu((i-1)*np*m + 2*np + 1 : (i-1)*np*m + 3*np, ...
        n*N + (i-1)*m*np - (i-1)*m + 1 : n*N + i*m*np - (i-1)*m) = blkdiag(Tcell_tauw{:});
    Tu((i-1)*np*m + 3*np + 1 : (i-1)*np*m + 4*np, ...
        n*N + (i-1)*m*np - (i-1)*m + 1 : n*N + i*m*np - (i-1)*m) = blkdiag(Tcell_tau{:});
    
end

auxdata.trans.T = T;
auxdata.trans.T_cost = T;
auxdata.trans.Tu = Tu;

end

