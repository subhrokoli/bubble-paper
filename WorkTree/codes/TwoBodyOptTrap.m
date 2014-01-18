% Code to solve the two-body optical trap Stokes flow problem
% Collaboration between IISER-K & IMSc
% 
% The first particle is trapped optically at the origin 
% and forced. 
% The second particle is free, and interacts with the first particle purely 
% hydrodynamically. The backreaction affects the first particle.
% 
% function [t X] = TwoBodyOptTrap(eta, k, a, d)
% 
% Inputs:
% 
% eta   :   fluid viscosity
% k     :   trap (harmonic) strength
% a     :   sphere radius
% d     :   separation of particles
% HI    :   switch for hydrodynamics. HI = 'on' or 'off'.
% 
% Outputs :
% 
% t     :   timeseries
% X     :   combined position vector of the particles. X(1:3) are the
%           position of the trapped particle, while X(4:6) are those of 
%           the free particle


% ======================================================================= %
% INTEGRATOR

function [t X] = TwoBodyOptTrap(eta, k, a, d, HI)

% Clock start
tic;

% compute trap force factor
trapFactor = k / (6*pi*eta*a);

% initial position of trapped particle
X10 = 1; Y10 = 0; Z10 = 0;         

% initial position of free particle
X20 = X10 + d; Y20 = Y10; Z20 = Z10;    

X0 = [X10 Y10 Z10 X20 Y20 Z20];

% timesteps and range
dt = 0.01;
t0 = 0.0; tf = 200.0;
tspan = t0:dt:tf; %lt = length(tspan);

% INTEGRATE
[t X] = ode15s(@(t,y)TwoBodyOptTrapDriver(t, x, trapFactor, a, HI), tspan, X0);

% Clock end
toc;

% EoF
end

