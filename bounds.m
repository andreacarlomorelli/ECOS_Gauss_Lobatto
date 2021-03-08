function auxdata = bounds(auxdata)
%
% B
% 
% FUNCTION DESCRIPTION: this function defines the boundaries of the
% variables
% 
% INPUTS/OUTPUTS:
% 
% auxdata     Structure containing auxiliary data 
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 01/03/2021
%

% Upper and lower bounds for state and control 
auxdata.bounds.lb = [0.1 -10 -10 1e-5 -10 0 log(0.1) -10 -10 -10 0];
auxdata.bounds.ub = [10 10 10 10 10 1000 0 10 10 10 10];

end