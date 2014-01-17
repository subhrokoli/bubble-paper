%% IISERK two particle data

%% Extract data
clear; clc
M = csvread('hyd_interact_4_1dot0.csv');

%% Plot

pl.hfig = figure; pp.OP = get( pl.hfig, 'outerposition' );
pp.scFac = 1.5;
set( pl.hfig, 'outerposition', [1 1 pp.scFac pp.scFac].*pp.OP );
pp.fs = 18;     % fontsize

plot(M(:, 1), M(:, 3), '.')
hl(1) = ylabel('Displacement of test bead (in Volt)');
hl(2) = xlabel('Time (sec)');
hl(3) = title('Laser amplitude : 1.0 V');

set(hl, 'fontsize', pp.fs)
set(gca, 'fontsize', pp.fs)
savePlot('hyd_interact_4_1dot0')