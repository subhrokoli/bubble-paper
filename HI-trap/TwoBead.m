%%
a = 1.5*1e-6;
y = 2.3*a;
eta = 0.001;

fc1 = 65; 
fc2 = 55;

k1 = 12*pi*pi*eta*a*fc1;
k2 = 12*pi*pi*eta*a*fc2;
f0 = 0.5*a;

% Parameters in terms of trap stiffness and other parameters
twoLambda = (k1+k2)/(6*pi*eta*a);
wo2 = k1*k2/(36*pi*pi*eta*eta*a*a) - k1*k2*(2 -4*a*a/(3*y*y))*(2 -4*a*a/(3*y*y))/(64*pi*pi*eta*eta*y*y);


% Amplitude of second bead 
A = f0*(2 -4*a*a/(3*y*y))/(8*pi*eta*y);
w = linspace(0.1,4000,256);
f_w =  A*w./sqrt((w.*w - wo2).*(w.*w - wo2) + twoLambda*twoLambda*w.*w);

aa = (wo2 - w.*w);
bb = twoLambda*w;
cc = wo2;
dd = (twoLambda - k2/(6*pi*eta*a))*w;

% Phase of the beads
delta1 = atan((aa.*dd - bb.*cc)./(aa.*cc + bb.*dd));
% Takes care of -90 to +90 jump in the phase.
% delta1 = pi*(delta1<0) +delta1; 

delta2 = atan(aa./bb);
% delta1 = atan(-w./(k1/(6*pi*eta*a))); Only first bead
delPhi = (delta2-delta1)*180/pi;

% Figure parameter
fig = figure(12);clf;
% set(gcf,'Renderer','zbuffer')

set(gcf,'OuterPosition', [993    55   930   993]);
% set(gca,'outerposition',[0.009 0.0 0.93 1]);
set(gca,'outerposition',[0.0 0.0 1 1]);
set(gca,'position',[0.13 0.13 0.81 0.81]);
set(gca,'fontsize',18);
set(0,'DefaultFigureColor','w',...
    'DefaultAxesColor',[1 1 1],...
    'DefaultAxesXColor',0.1*[1 1 1],...
    'DefaultAxesYColor',0.1*[1 1 1],...
    'DefaultAxesZColor',0.1*[1 1 1],...
    'DefaultTextColor','k',...
    'DefaultLineColor','k')

plot(w/(2*pi), -delta1*180/pi,'-','LineWidth',2,'Color',[154 0 0]./256,'MarkerSize',10); hold on;
plot(w/(2*pi), -delta2*180/pi,'-','LineWidth',2,'Color',[0 134 212]./256,'MarkerSize',10);
plot(w/(2*pi), -delPhi, '-','LineWidth',2,'Color',[98 181 3]./256,'MarkerSize',10);grid on;
legend({'$\delta_1$','$\delta_2$','$\Delta \delta$'},'FontSize',16,'interpreter','latex','location','southeast');
grid on;
figure(2);
plot(w/(2*pi),f_w,'LineWidth',2,'Color',[0 134 212]./256);line([77,77],[0,16e-8],'LineWidth',1.2,'Color',[154 0 0]./256);
% line([77,77],[-100,150],'LineWidth',1,'Color',[96 0 148]./256);
grid on;
% hl(1) = xlabel('Driving Frequency (Hz)');
% hl(2) = ylabel('Amplitude of second particle (in micron)');
set(hl,'fontsize',18);
axis square;box on;