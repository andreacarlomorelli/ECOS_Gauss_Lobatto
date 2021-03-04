function auxdata = Earth_Venus_parameters
%
% Earth_Venus_parameters
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
auxdata.engine.Tmax = 0.33; % Maximum thrust
auxdata.engine.Isp = 3800; % Specific Impulse
auxdata.sc.m0 = 1500; % Initial mass
auxdata.engine.c = auxdata.engine.Tmax*auxdata.units.R0/(auxdata.sc.m0*auxdata.units.V0^2); % Constant c
auxdata.engine.ve = auxdata.engine.Isp*g0; % Exhaust velocity
auxdata.t0 = 0; % Adimensionalized initial time
auxdata.tf = 1000*24*3600/(auxdata.units.R0/auxdata.units.V0); % Adimesionalized final time

auxdata.rr_0 = [0.9708, 0.2376, 0]'; % Initial position
auxdata.vv_0 = [-0.2545, 0.9687, 0]'; % Initial velocity

auxdata.rr_f = [-0.3277, 0.6389, 0.0277]'; % Final position
auxdata.vv_f = [-1.0509, -0.5436, 0.0532]'; % Final velocity

auxdata.eps_x0 = 1e-8;
auxdata.eps_xf = 1e-8;
end
