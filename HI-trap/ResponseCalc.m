%%
close all;
figure(12); clf;
eta = 1/6;
a   = 1;
y = 5*a;

mu = 1/(6*pi*eta*a);
mu12 = (2 -4*a*a/(3*y*y))/(8*pi*eta*y);

k1 = 20;
k2 = 10;

DetA = mu^2*(k1*k2) - mu12^2*k1*k2;
TrA  = mu*(k1 + k2);

tau1 = 1/(mu*k1);
tau2 = 1/(mu*k2);

w = linspace(-40,40,1024);

%%
chi_p11 = -(mu*k2*(DetA-w.^2) + w.^2*TrA)./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p12 =  (mu12*k1*(DetA-w.^2))         ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p21 =  (mu12*k2*(DetA-w.^2))         ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p22 = -(mu*k1*(DetA-w.^2) + w.^2*TrA)./((DetA-w.^2).^2 + w.^2*TrA.^2);

chi_2p11 = -w.*( DetA-w.^2 -mu*k2*TrA )./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p12 = -w.*( mu12*k1*TrA)          ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p21 = -w.*( mu12*k2*TrA)          ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p22 = -w.*( DetA-w.^2 -mu*k1*TrA )./((DetA-w.^2).^2 + w.^2*TrA.^2);


wp = w*tau2;
xp = chi_p12;
yp = chi_2p12;

plot(wp,xp,'*','color', [0.0 0.5 0.9]);hold on;
plot(wp,yp,'*','color', [0.8 0.1 0.1]);
xlabel('$\omega \tau$','interpreter','latex');
ylabel('$\chi$','interpreter','latex');
axis tight;
legend({'$\chi^{\prime}$','$\chi^{\prime\prime}$'},'interpreter','latex');
figure(13);
plot(xp,yp,'-*','color', [0.4 0.7 0.4]);
axis tight;
figure(23);hold on;
plot(wp, xp.*xp + yp.*yp,'*','color', [0.2 0.6 0.5]);
% plot(wp,atan(yp./xp),'*','color', [0.2 0.6 0.5]);
axis tight;

%%

chi_2p11 = -w.*(mu*k2*(DetA-w.^2) + w.^2*TrA)./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p12 =  w.*(mu12*k1*(DetA-w.^2))         ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p21 =  w.*(mu12*k2*(DetA-w.^2))         ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p22 = -w.*(mu*k1*(DetA-w.^2) + w.^2*TrA)./((DetA-w.^2).^2 + w.^2*TrA.^2);

chi_p11 = w.*w.*( DetA-w.^2 -mu*k2*TrA )./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p12 = w.*w.*( mu12*k1*TrA)          ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p21 = w.*w.*( mu12*k2*TrA)          ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p22 = w.*w.*( DetA-w.^2 -mu*k1*TrA )./((DetA-w.^2).^2 + w.^2*TrA.^2);


wp = w*tau2;
xp = chi_p12;
yp = chi_2p12;

plot(wp,xp,'*','color', [0.0 0.5 0.9]);hold on;
plot(wp,yp,'*','color', [0.8 0.1 0.1]);
xlabel('$\omega \tau$','interpreter','latex');
ylabel('$\chi$','interpreter','latex');
axis tight;
legend({'$\chi^{\prime}$','$\chi^{\prime\prime}$'},'interpreter','latex');
figure(13);
plot(xp,yp,'-*','color', [0.4 0.7 0.4]);
axis tight;
figure(23);hold on;
plot(wp, xp.*xp + yp.*yp,'*','color', [0.2 0.6 0.5]);
% plot(wp,atan(yp./xp),'*','color', [0.2 0.6 0.5]);
axis tight;
