function auxdata = Earth_Dionysus_parameters
%
% Earth_Dionysus_parameters
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
auxdata.engine.Tmax = 0.32; % Maximum thrust
auxdata.engine.Isp = 3000; % Specific Impulse
auxdata.sc.m0 = 4000; % Initial mass
auxdata.engine.c = auxdata.engine.Tmax*auxdata.units.R0/(auxdata.sc.m0*auxdata.units.V0^2); % Constant c
auxdata.engine.ve = auxdata.engine.Isp*g0; % Exhaust velocity
auxdata.t0 = 0; % Adimensionalized initial time
auxdata.tf = 3534*24*3600/(auxdata.units.R0/auxdata.units.V0); % Adimesionalized final time

auxdata.rr_0 = [-0.0243, 0.9833, -1.5117e-5]'; % Initial position
auxdata.vv_0 = [-1.0161, -0.0285, 1.6955e-6]'; % Initial velocity

auxdata.rr_f = [-2.0406, 2.0518, 0.5543]'; % Final position
auxdata.vv_f = [-0.1423, -0.4511, 0.0189]'; % Final velocity

auxdata.eps_x0 = 1e-8;
auxdata.eps_xf = 1e-8;

end
