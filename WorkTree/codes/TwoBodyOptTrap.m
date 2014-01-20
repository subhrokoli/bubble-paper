% ======================================================================= %
% MAIN CODE
% Code to solve the two-body optical trap Stokes flow problem
% Collaboration between IISER-K & IMSc
% 
% The first particle is trapped optically at the origin using a harmonic
% potential. It is sinusoidally forced.
% The second particle is free, and interacts with the first particle purely 
% hydrodynamically. The backreaction affects the first particle.
% 
% function [t R] = TwoBodyOptTrap(eta, k, F0, omega, a, R0)
% 
% Inputs:
% 
% eta   :   fluid viscosity
% k     :   trap (harmonic) strength
% F0    :   magnitude of sinusoidal driving force
% omega :   frequency of sinusoidal driving force
% a     :   sphere radius
% R0    :   initial position of particles (6 X 1 vector)
% 
% Outputs :
% 
% t     :   timeseries
% R     :   combined position vector of the particles. R(1:3) are the
%           position of the trapped particle, while R(4:6) are those of 
%           the free particle
% 
% (c) Somdeb Ghose, IMSc, Chennai, India (2014)
% ======================================================================= %

function [t R] = TwoBodyOptTrap(eta, k, F0, omega, a, R0)

% Clock start
tic;

% compute mobility factor
mob_fac = 1 / (6*pi*eta*a);

% % set initial position of particles
% % first three are coordinates of trapped particle
% % second three are coordinates of free particle
% R0 = [1, 0, 0, d, 0, 0];

% initialize stresslet and save
S = zeros(18, 1); save('S.mat', 'S');

% timesteps and range
dt = 0.01;
t0 = 0.0; tf = 200.0;
tspan = t0:dt:tf; %lt = length(tspan);

% INTEGRATE
[t R] = ode15s(@(t,r)TwoBodyDriver(t, r, a, mob_fac, k, F0, omega), tspan, R0);

% plot
h.plot = plot(t, R);
grid on; 

% Clock end
toc;

% EoF
end