function auxdata = Earth_Mars_parameters 
%
% Earth_Mars_parameters
% 
% FUNCTION DESCRIPTION: this function defines the parameters for the Earth
% to Dionysus stransfer
% 
% INPUTS/OUTPUTS:
% 
% auxdata     Structure containing auxiliary data 
%
% AUTHOR: Andrea Carlo Morelli
% 
% DATE: 01/03/2021
%

% Physical constants
g0 = 9.80665;
auxdata.units.R0 = 1.49597870e11; % Departure radius: Earth
auxdata.gravitational.mu = 1.3271244e20; % Gravitational constant of the Sun
auxdata.units.V0 = sqrt(auxdata.gravitational.mu/auxdata.units.R0); % Velocity of the Earth

% Transfer parameters
auxdata.engine.Tmax = 0.5; % Maximum thrust
auxdata.engine.Isp = 2000; % Specific Impulse
auxdata.sc.m0 = 1000; % Initial mass
auxdata.engine.c = auxdata.engine.Tmax*auxdata.units.R0/(auxdata.sc.m0*auxdata.units.V0^2); % Constant c
auxdata.engine.ve = auxdata.engine.Isp*g0; % Exhaust velocity
auxdata.t0 = 0; % Adimensionalized initial time
auxdata.tf = 348.795*24*3600/(auxdata.units.R0/auxdata.units.V0); % Adimesionalized final time

auxdata.rr_0 = [-140699693, -51614428, 980]'*1000/auxdata.units.R0; % Initial position
auxdata.vv_0 = [9.774596, -28.07828, 4.337725e-4]'*1000/auxdata.units.V0; % Initial velocity

auxdata.rr_f = [-172682023, 176959469, 7948912]'*1000/auxdata.units.R0; % Final position
auxdata.vv_f = [-16.427384, -14.860506, 9.21486e-2]'*1000/auxdata.units.V0; % Final velocity

auxdata.eps_x0 = 1e-8;
auxdata.eps_xf = 1e-8;

end
