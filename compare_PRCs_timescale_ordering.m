% ==================================================================
%
%                    compare_PRCs_timescale_ordering.m
%                      ------- 
%  Compares the PRCs for the neuromechanical module
%  in different timescale orderings

%  Makes several figures:
%       PRCs in three parameter regimes
%      
% ==================================================================
addpath('./src');
clear
load('colorblind_colormap.mat');

% ------ tau_b > tau_m ----------------
tic
% -- NM MODEL PARAMETERS --
mu = 1.3e-7; %my larger mu so that mu/kb = 0.5
kb = 2.6e-7;
t_b=mu/kb 

t_n = 1e-2; t_m = 1.5e-1; %timescales for length, neural, and muscule activity

c_ma = 5; c_prop = 1;  %musc. activity feedback strength, prop feedback strength
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 2; %chain of 2 units

little_gamma = [1, 10, logspace(2,3,4), 10^4, 10^4.4472];
eps_prop = 0.05;
eps_gap = 0.015;

% time step size
dt=1e-3;

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_b, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_b, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(1); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, \tau_b > \tau_m');
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);

% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig4 = figure(4); clf;
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) > 0);
unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25); 

%figure 5: compare curvature PRC, coupling current (K'_LC), and mech G-fn
fig5 = figure(5); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('\tau_b > \tau_m'); 
subplot(3,1,2);
%approximate K' from limit cycle with finite differences
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);
plot(0:dt/period:1-dt/period, Kdv,'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('d\kappa/dt'); set(gca,'FontSize',30);
% %plot Gfn
% subplot(3,1,3)
% plot(theta,g1,'r', 'linewidth',4); hold on;
% plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
%     0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
%     'Markersize',15);
% plot(theta(intersect(pl_state_inds, stab_state_inds)),...
%     0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
%     'Markersize',20, 'MarkerFaceColor', 'k');
% plot([0,1],[0,0],'k:','linewidth',2); hold off;
% xlabel('\theta', 'FontSize', 25); ylabel('G_{m}');
% set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);
subplot(3,1,3)
hm1 = plot(theta,-H_m(-theta),'-r', 'linewidth',4); hold on;
hm2 = plot(theta,-H_m(theta),'-.r', 'linewidth',4); hold on;
g = plot(theta,g1,':r', 'linewidth',4); hold on;
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('\theta', 'FontSize', 25); ylabel('H_{m}');
% legend([hm1 hm2], 'H_m(-\phi)', 'H_m(\phi)')
legend([hm1 hm2, g], {'H_m(-\phi)', 'H_m(\phi)', 'G_m(\phi)'})
set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);


%show with Sinusoidal coupling current
% ---- III.  CALCULATE G-FUNCTION  ----
%replace X with cosine so that mech. current is sin
X(:,1) = 10*period/(2*pi).*cos(2*pi.*(dt:dt:period)./period);
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig11 = figure(11); clf;
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) > 0);
unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25); 

%figure 12: compare curvature PRC, coupling current (K'_LC), and mech G-fn
fig12 = figure(12); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('\tau_b > \tau_m'); 
%approximate K' from limit cycle with finite differences
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);
subplot(3,1,2);
plot(0:dt/period:1-dt/period, Kdv,'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('d\kappa/dt'); set(gca,'FontSize',30);
%plot H fns
subplot(3,1,3)
hm1 = plot(theta,-H_m(-theta),'-r', 'linewidth',4); hold on;
hm2 = plot(theta,-H_m(theta),'-.r', 'linewidth',4);
g = plot(theta,g1,':r', 'linewidth',4);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('\theta', 'FontSize', 25); ylabel('H_{m}');
legend([hm1 hm2 g], {'H_m(-\phi)', 'H_m(\phi)', 'G_m(\phi)'})
set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);


 
% ------ tau_b < tau_m ----------------
% ----------------------------------
t_b = 0.05; %t_b < t_m=0.15

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_b, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_b, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(2); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, \tau_b < \tau_m'); 
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-.2,.2]);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);


% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig6 = figure(6); clf;
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) > 0);
unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)

%figure 7: compare curvature PRC, coupling current (K'_LC), and mech G-fn
fig7 = figure(7); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('\tau_b < \tau_m'); 
subplot(3,1,2);
%approximate K' from limit cycle with finite differences
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);
plot(0:dt/period:1-dt/period, Kdv,'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('d\kappa/dt'); set(gca,'FontSize',30);
% %plot Gfn
% subplot(3,1,3)
% plot(theta,g1,'r', 'linewidth',4); hold on;
% plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
%     0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
%     'Markersize',15);
% plot(theta(intersect(pl_state_inds, stab_state_inds)),...
%     0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
%     'Markersize',20, 'MarkerFaceColor', 'k');
% plot([0,1],[0,0],'k:','linewidth',2); hold off;
% xlabel('\theta', 'FontSize', 25); ylabel('G_{m}');
% set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);
%plot Hfns
subplot(3,1,3)
hm1 = plot(theta,-H_m(-theta),'-r', 'linewidth',4); hold on;
hm2 = plot(theta,-H_m(theta),'-.r', 'linewidth',4); hold on;
g = plot(theta,g1,':r', 'linewidth',4); hold on;
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('\theta', 'FontSize', 25); ylabel('H_{m}');
% legend([hm1 hm2], 'H_m(-\phi)', 'H_m(\phi)')
legend([hm1 hm2, g], {'H_m(-\phi)', 'H_m(\phi)', 'G_m(\phi)'})
set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);

%show with Sinusoidal coupling current
% ---- III.  CALCULATE G-FUNCTION  ----
%replace X with cosine so that mech. current is sin
X(:,1) = 10*period/(2*pi).*cos(2*pi.*(dt:dt:period)./period);
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig13 = figure(13); clf;
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) > 0);
unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25); 

%figure 14: compare curvature PRC, coupling current (K'_LC), and mech G-fn
fig14 = figure(14); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('\tau_b < \tau_m'); 
%approximate K' from limit cycle with finite differences
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);
subplot(3,1,2);
plot(0:dt/period:1-dt/period, Kdv,'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('d\kappa/dt'); set(gca,'FontSize',30);
%plot H fns
subplot(3,1,3)
hm1 = plot(theta,-H_m(-theta),'-r', 'linewidth',4); hold on;
hm2 = plot(theta,-H_m(theta),'-.r', 'linewidth',4); hold on;
g = plot(theta,g1,':r', 'linewidth',4); hold on;
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('\theta', 'FontSize', 25); ylabel('H_{m}');
legend([hm1 hm2, g], {'H_m(-\phi)', 'H_m(\phi)', 'G_m(\phi)'})
set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);

return

% ------ tau_b \approx tau_m ----------------
% ----------------------------------
t_b = 0.14; %t_b \approx t_m=0.15

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_b, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_b, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(3); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, \tau_b \approx \tau_m'); 
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);


% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig8 = figure(8); clf;
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) > 0);
unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)

%figure 9: compare curvature PRC, coupling current (K'_LC), and mech G-fn
fig9 = figure(9); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('\tau_b \approx \tau_m'); 
subplot(3,1,2);
%approximate K' from limit cycle with finite differences
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);
plot(0:dt/period:1-dt/period, Kdv,'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('d\kappa/dt'); set(gca,'FontSize',30);
%plot Gfn
subplot(3,1,3)
plot(theta,g1,'r', 'linewidth',4); hold on;
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
xlabel('\theta', 'FontSize', 25); ylabel('G_{m}');
set(gca,'FontSize',30); set(gcf,'Position',[1    59   640   646]);


%Make H-
