%%
close all;
eta = 1/6;
a   = 1;
y = 2*a;

mu = 1/(6*pi*eta*a);
mu12 = (2 -4*a*a/(3*y*y))/(8*pi*eta*y);

k1 = 40;
k2 = 50;

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


wp = w*tau1;
xp = chi_p11;
yp = chi_2p11;


figure(12); clf;

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

chi_p11 =  (k2*DetA*(mu^2-mu12^2) + w.*w.*(k2*mu12^2 +k1*mu^2))./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p12 =  w.*w.*TrA*mu12                                      ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p21 =  w.*w.*TrA*mu12                                      ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p22 =  (k1*DetA*(mu^2-mu12^2) + w.*w.*(k1*mu12^2 +k2*mu^2))./((DetA-w.^2).^2 + w.^2*TrA.^2);

chi_2p11 = w.*( mu*k2*(mu^2-mu12^2) + mu*w.*w)./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p12 = -w.*mu12.*(DetA - w.*w)             ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p21 = -w.*mu12.*(DetA - w.*w)             ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p22 = w.*( mu*k2*(mu^2-mu12^2) + mu*w.*w)./((DetA-w.^2).^2 + w.^2*TrA.^2);


wp = w*tau1;
xp = chi_p21;
yp = chi_2p21;


figure(117); clf;
set(gcf,'outerposition', [400 400 1400 550]);
% set(gca,'outerposition',[0.009 0.0 0.93 1]);
% set(gca,'position',[0.1 0.1 0.85 0.85]);
set(0,'DefaultFigureColor','w',...
    'DefaultAxesColor',[1 1 1],...
    'DefaultAxesXColor',0.1*[1 1 1],...
    'DefaultAxesYColor',0.1*[1 1 1],...
    'DefaultAxesZColor',0.1*[1 1 1],...
    'DefaultTextColor','k',...
    'DefaultLineColor','k');
h1 = subplot(1,3,1);
set(h1,'position',[0.06 0.12 0.26 0.82]);
plot(wp,xp,'-','color', [0.0 0.5 0.9],'linewidth', 2);hold on;
plot(wp,yp,'-','color', [0.8 0.1 0.1],'linewidth', 2);
xlabel('\boldmath{$\Omega \tau_1$}','interpreter','latex','fontsize',16);
ylabel('\boldmath{$\chi_{21}$}','interpreter','latex','fontsize',16);
axis tight;
axis square;
grid on;
set(gca,'fontsize', 16);
legend({'\boldmath {$\chi_{21}^{\prime}$}','\boldmath{$\chi_{21}^{\prime\prime}$}'},'interpreter','latex');
% figure(117);

h2 = subplot(1,3,2);
set(h2,'position',[0.39 0.12 0.26 0.82]);
plot(xp(512:end),yp(512:end),'-','color', [0.4 0.7 0.4],'linewidth', 2);
axis tight;
axis square;
xlabel('\boldmath{$\chi_{21}^{\prime}$}','interpreter','latex','fontsize',16);
ylabel('\boldmath{$\chi_{21}^{\prime \prime}$}','interpreter','latex','fontsize',16);
set(gca,'fontsize', 16);
grid on;
arrow([0.00038 -0.0012],[0.001 -0.0022],'linewidth',1.6, 'baseangle', 15)
text(0.0007,-0.0015,'\boldmath $\Omega \tau_1$','interpreter','latex','fontsize', 16);

h3 = subplot(1,3,3);
set(h3,'position',[0.72 0.12 0.26 0.82]);
% figure(119);
hold on;
plot(wp, sqrt(xp.*xp + yp.*yp),'-','color', [0.2 0.6 0.5],'linewidth', 2);
axis tight;box on;
ylabel('\boldmath{$|\chi_{21}|$}','interpreter','latex','fontsize',16);
xlabel('\boldmath{$\Omega \tau_1$}','interpreter','latex','fontsize',16);
set(gca,'fontsize', 16);
grid on;

axis square;
% figure(123);
% plot(wp,atan(-yp./xp),'*','color', [0.2 0.6 0.5]);

%% CHI11


chi_p11 =  (k2*DetA*(mu^2-mu12^2) + w.*w.*(k2*mu12^2 +k1*mu^2))./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p12 =  w.*w.*TrA*mu12                                      ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p21 =  w.*w.*TrA*mu12                                      ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_p22 =  (k1*DetA*(mu^2-mu12^2) + w.*w.*(k1*mu12^2 +k2*mu^2))./((DetA-w.^2).^2 + w.^2*TrA.^2);

chi_2p11 = w.*( mu*k2*(mu^2-mu12^2) + mu*w.*w)./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p12 = -w.*mu12.*(DetA - w.*w)             ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p21 = -w.*mu12.*(DetA - w.*w)             ./((DetA-w.^2).^2 + w.^2*TrA.^2);
chi_2p22 = w.*( mu*k2*(mu^2-mu12^2) + mu*w.*w)./((DetA-w.^2).^2 + w.^2*TrA.^2);


wp = w*tau1;
xp = chi_p22;
yp = chi_2p22;


figure(117); clf;
set(gcf,'outerposition', [400 400 1400 550]);
% set(gca,'outerposition',[0.009 0.0 0.93 1]);
% set(gca,'position',[0.1 0.1 0.85 0.85]);
set(0,'DefaultFigureColor','w',...
    'DefaultAxesColor',[1 1 1],...
    'DefaultAxesXColor',0.1*[1 1 1],...
    'DefaultAxesYColor',0.1*[1 1 1],...
    'DefaultAxesZColor',0.1*[1 1 1],...
    'DefaultTextColor','k',...
    'DefaultLineColor','k');
h1 = subplot(1,3,1);
set(h1,'position',[0.06 0.12 0.26 0.82]);
plot(wp,xp,'-','color', [0.0 0.5 0.9],'linewidth', 2);hold on;
plot(wp,yp,'-','color', [0.8 0.1 0.1],'linewidth', 2);
xlabel('\boldmath{$\Omega \tau_1$}','interpreter','latex','fontsize',16);
ylabel('\boldmath{$\chi_{11}$}','interpreter','latex','fontsize',16);
axis tight;
axis square;
grid on;
set(gca,'fontsize', 16);
legend({'\boldmath {$\chi_{11}^{\prime}$}','\boldmath{$\chi_{11}^{\prime\prime}$}'},'interpreter','latex');
% figure(117);

h2 = subplot(1,3,2);
set(h2,'position',[0.39 0.12 0.26 0.82]);
plot(xp(512:end),yp(512:end),'-','color', [0.4 0.7 0.4],'linewidth', 2);
axis tight;
axis square;
xlabel('\boldmath{$\chi_{11}^{\prime}$}','interpreter','latex','fontsize',16);
ylabel('\boldmath{$\chi_{11}^{\prime \prime}$}','interpreter','latex','fontsize',16);
set(gca,'fontsize', 16);
grid on;
arrow([0.0175 0.00045],[0.0154 0.0013],'linewidth',1.6, 'baseangle', 15)
text(0.015,0.0006,'\boldmath $\Omega \tau_1$','interpreter','latex','fontsize', 16);

h3 = subplot(1,3,3);
set(h3,'position',[0.72 0.12 0.26 0.82]);
% figure(119);
hold on;
plot(wp, sqrt(xp.*xp + yp.*yp),'-','color', [0.2 0.6 0.5],'linewidth', 2);
axis tight;box on;
ylabel('\boldmath{$|\chi_{11}|$}','interpreter','latex','fontsize',16);
xlabel('\boldmath{$\Omega \tau_1$}','interpreter','latex','fontsize',16);
set(gca,'fontsize', 16);
grid on;

axis square;
% figure(123);
% plot(wp,atan(-yp./xp),'*','color', [0.2 0.6 0.5]);
