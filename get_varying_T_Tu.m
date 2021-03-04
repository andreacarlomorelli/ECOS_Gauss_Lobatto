function [paraECOS, paraGL, paraSCP, auxdata] = get_varying_T_Tu(x_old, paraECOS, paraGL, paraSCP, auxdata)
%
% get_varying_T_Tu
% 
% FUNCTION DESCRIPTION: this function computes the varying parts of the
% matrices T and Tu
% 
% INPUTS:
%
% x_old          Reference solution                          
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
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 25/02/2021
%

% Sates and controls
n = 7;
m = 4;

% GL parameters
np = paraGL.np;
Ni = paraGL.Ni;
N = paraGL.N;

% SCP parameters
h = paraSCP.h;

% Useful parameters
c = auxdata.engine.c;
ve = auxdata.engine.ve;
V0 = auxdata.units.V0;

% Matrix A for all the nodes
As = zeros(n,n,N);
for i = 1 : N
    As(:,:,i) = A(x_old(i,1:7));
end
paraECOS.As = As;

% Matrix B
Bs = B(c, ve, V0);

% T matrices
T = auxdata.trans.T;

for i = 1 : Ni
    for k = 1 : n
        Abig = zeros(np,n*np);
        Bbig = zeros(np,m*np);
        for j = 1 : np
            Abig(j, (j-1)*n + 1 : j*n) = h*0.5*As(k,:,(i - 1) * (np - 1) + j);
            Bbig(j, (j-1)*m + 1 : j*m) = h*0.5*Bs(k,:);  
        end
        T((i-1)*2*np*n + 2*(k-1)*np + np + 1: (i-1)*2*np*n + ...
            2*(k-1)*np + 2*np, (i-1)*n*np - (i-1)*n + 1 : i*n*np - (i-1)*n) = Abig;
        T((i-1)*2*np*n + 2*(k-1)*np + np + 1: (i-1)*2*np*n + ...
            2*(k-1)*np + 2*np, N*n + (i-1)*m*np + 1 - (i-1)*m ...
            : N*n + i*m*np - (i-1)*m) = Bbig;
    end
end

auxdata.trans.T = T;

end
