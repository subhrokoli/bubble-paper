clear all;close all;clc
filename = 'water_calibrated_brownian.xlsx';
sheet = 1;
xlRange = 'B2:B20001';
x = xlsread(filename,sheet,xlRange);
N = length(x);
x_n1 = x(2:N)'*x(2:N);
x_n = x(1:N-1)'*x(1:N-1);
x_pairs = x(1:N-1)'*x(2:N);
x_1 = (x(1))^2;

dt = 0.5e-03;

LAMBDA = linspace(302,306,100);
d = linspace(0.3,0.5,100);

[lambda,D] = meshgrid(LAMBDA,d);

a = 1 - exp(-2*lambda*dt);
b = exp(-lambda*dt);
c = exp(-2*lambda*dt);

lnP = (((N-1)/2)*log(lambda./D)) - (((N-1)/2)*log(2*pi*a)) - ((lambda./D).*(1./(2*a))*(x_n1 + x_n.*c -2*x_pairs.*b)) + (0.5*log(lambda./(2*pi*D))) - (lambda*x_1./(2*D));
surf(lambda,D,lnP)