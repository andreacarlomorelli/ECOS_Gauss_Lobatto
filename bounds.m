function auxdata = bounds(auxdata)
%
% bounds
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
auxdata.bounds.lb = [-10 -10 -10 -10 -10 -10 log(0.1) -10 -10 -10 0];
auxdata.bounds.ub = [10 10 10 10 10 10 0 10 10 10 10];

end