% ==================================================================
%
%                     phase_response_properties.m
%
%  This program plots the phase response properties 
%   and G-functions (why not?)
%
% ==================================================================

addpath('./src');
clear
load('colorblind_colormap.mat');

%--------------------------------------------------------------------
%Part 1 - Parameters from Paper (ch2 and 3 of dissertation) (tristable
%regime)

% -- Fixed NM MODEL PARAMETERS --
mu = 1.3e-7;
kb = 2.6e-7;
t_f=mu/kb;
t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
c_m = 10; c_p = 1;
a = 1; I = 0; %neural voltage model param, AVB input bias current
sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
dt=1e-3; % time step size


% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_m, c_p, t_f, t_n, t_m, a, I, sigma );

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_m, c_p, t_f, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(1); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, Tristable Regime'); 
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);



% ---- III.  CALCULATE G-FUNCTION  ----
% -- Mechanical coupling -computing weight matrix B---
gridsz=1; dim =2;
delX = 1/(gridsz*dim);
e = ones(gridsz*dim,1);
%2nd difference operator
D2 = spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%4th difference operator
D4 = D2*D2';
D4(1,1) = 7; D4(end,end) = 7; %moment-free, force-free BCs
D4 = (1/delX^4).*D4;
%mechanical coupling weight matrix B
B = inv(full(D4));

%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g_gap = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig = figure(10); clf;
% subplot(2,1,1); 
mech1 = plot(theta,g1,'-.r', 'linewidth',4); hold on;
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
plot([0,1],[0,0],'k:','linewidth',2);
% title('\tau_b = 10/3 \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25);
ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)

%G fn - gap junctional coupling
fig2 = figure(12); clf;
gap1 = plot(theta,g_gap/max(g_gap),'-.b', 'linewidth',4); hold on; 
%find zeros & mark stability
pl_state_inds = find(abs(g_gap)/max(g_gap)<=1e-3);
stab_state_inds = find(g_gap-circshift(g_gap,-1) >= 0);
unstab_state_inds = find(g_gap-circshift(g_gap,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2);
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_g(\phi)'); set(gca,'FontSize',30)

%G fn - proprioceptive coupling
fig3 = figure(13); load('colorblind_colormap.mat');
prop1 = plot(theta,h1p/max(h1p),'-.', 'linewidth',4, 'Color', colorblind(10,:,:)); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(h1p)/max(h1p)<=1e-3);
stab_state_inds = find(h1p-circshift(h1p,-1) >= 0);
unstab_state_inds = find(h1p-circshift(h1p,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); 
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_p(\phi)'); set(gca,'FontSize',30)

%save Z for comparisons later
Z2 = Z;
%--------------------------------------------------------------------
%Part 2 - Bistable regime

% -- Fixed NM MODEL PARAMETERS --
t_f = 0.5;
t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
c_m = 10; c_p = 10;
a = 1; I = 0.5; %neural voltage model param, AVB input bias current
sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
dt=1e-3; % time step size

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_m, c_p, t_f, t_n, t_m, a, I, sigma);

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_m, c_p, t_f, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(2); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, Bistable Regime'); 
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.1,0.1]);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.05,0.05]);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);

% ---- III.  CALCULATE G-FUNCTION  ----
% -- Mechanical coupling -computing weight matrix B---
gridsz=1; dim =2;
delX = 1/(gridsz*dim);
e = ones(gridsz*dim,1);
%2nd difference operator
D2 = spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%4th difference operator
D4 = D2*D2';
D4(1,1) = 7; D4(end,end) = 7; %moment-free, force-free BCs
D4 = (1/delX^4).*D4;
%mechanical coupling weight matrix B
B = inv(full(D4));

%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g_gap = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

% fig = figure(10); 
% % subplot(2,1,1); 
% m2 = plot(theta,g1,'r', 'linewidth',4); hold on;
% %find zeros & mark stability
% pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
% stab_state_inds = find(g1-circshift(g1,-1) > 0);
% unstab_state_inds = find(g1-circshift(g1,-1) <= 0);
% plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
%     0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
%     'Markersize',15);
% plot(theta(intersect(pl_state_inds, stab_state_inds)),...
%     0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
%     'Markersize',20, 'MarkerFaceColor', 'k');
% plot([0,1],[0,0],'k:','linewidth',2);
% % title('\tau_b = 10/3 \tau_m', 'FontSize',25);
% xlabel('phase difference \phi', 'FontSize', 25);
% ylabel('G_{m}(\phi)','FontSize', 25);
% set(gca,'FontSize',25)
% legend([m1 m2],'tristable','bistable')
% 
% %G fn - gap junctional coupling
% fig2 = figure(12);
% gap2 = plot(theta,g_gap/max(g_gap),'b', 'linewidth',4); hold on; 
% %find zeros & mark stability
% pl_state_inds = find(abs(g_gap)/max(g_gap)<=1e-3);
% stab_state_inds = find(g_gap-circshift(g_gap,-1) >= 0);
% unstab_state_inds = find(g_gap-circshift(g_gap,-1) < 0);
% plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
%     0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
%     'Markersize',15);
% plot(theta(intersect(pl_state_inds, stab_state_inds)),...
%     0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
%     'Markersize',20, 'MarkerFaceColor', 'k');
% plot([0,1],[0,0],'k:','linewidth',2);
% xlabel('phase \phi');  ylim([-1 1]); ylabel('G_g(\phi)'); set(gca,'FontSize',30)
% legend([gap1 gap2],'tristable','bistable')
% 
% %G fn - proprioceptive coupling
% fig3 = figure(13);  load('colorblind_colormap.mat');
% prop2 = plot(theta,h1p/max(h1p),'-', 'linewidth',4, 'Color', colorblind(10,:,:)); hold on;
% %find zeros & mark stability
% pl_state_inds = find(abs(h1p)/max(h1p)<=1e-3);
% stab_state_inds = find(h1p-circshift(h1p,-1) >= 0);
% unstab_state_inds = find(h1p-circshift(h1p,-1) < 0);
% plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
%     0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
%     'Markersize',15);
% plot(theta(intersect(pl_state_inds, stab_state_inds)),...
%     0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
%     'Markersize',20, 'MarkerFaceColor', 'k');
% plot([0,1],[0,0],'k:','linewidth',2); 
% xlabel('phase \phi');  ylim([-1 1]); ylabel('G_p(\phi)'); set(gca,'FontSize',30)
% % legend([prop1 prop2],'tristable','bistable')


%--------------------------------------------------------------------
%Part 3 - Effect of Cm

% -- Fixed NM MODEL PARAMETERS --
mu = 1.3e-7;
kb = 2.6e-7;
t_f=mu/kb;
t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
c_m = 20; c_p = 1;
a = 1; I = 0; %neural voltage model param, AVB input bias current
sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
dt=1e-3; % time step size


% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_m, c_p, t_f, t_n, t_m, a, I, sigma );

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_m, c_p, t_f, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(3); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, Increased C_m'); 
ylim([-0.2,0.2]);
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.1,0.1]);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.05,0.05]);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);


% ---- III.  CALCULATE G-FUNCTION  ----
% -- Mechanical coupling -computing weight matrix B---
gridsz=1; dim =2;
delX = 1/(gridsz*dim);
e = ones(gridsz*dim,1);
%2nd difference operator
D2 = spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%4th difference operator
D4 = D2*D2';
D4(1,1) = 7; D4(end,end) = 7; %moment-free, force-free BCs
D4 = (1/delX^4).*D4;
%mechanical coupling weight matrix B
B = inv(full(D4));

%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g_gap = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig = figure(10);
% subplot(2,1,1); 
mech2 = plot(theta,g1,'-r', 'linewidth',4); hold on;
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
plot([0,1],[0,0],'k:','linewidth',2);
% title('\tau_b = 10/3 \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25);
ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)

%G fn - gap junctional coupling
fig2 = figure(12);
gap2 = plot(theta,g_gap/max(g_gap),'-b', 'linewidth',4); hold on; 
%find zeros & mark stability
pl_state_inds = find(abs(g_gap)/max(g_gap)<=1e-3);
stab_state_inds = find(g_gap-circshift(g_gap,-1) >= 0);
unstab_state_inds = find(g_gap-circshift(g_gap,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2);
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_g(\phi)'); set(gca,'FontSize',30)

%G fn - proprioceptive coupling
fig3 = figure(13);
prop2 = plot(theta,h1p/max(h1p),'-', 'linewidth',4, 'Color', colorblind(10,:,:)); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(h1p)/max(h1p)<=1e-3);
stab_state_inds = find(h1p-circshift(h1p,-1) >= 0);
unstab_state_inds = find(h1p-circshift(h1p,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); 
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_p(\phi)'); set(gca,'FontSize',30)



%--------------------------------------------------------------------
%Part 3 - Effect of Cp

% -- Fixed NM MODEL PARAMETERS --
mu = 1.3e-7;
kb = 2.6e-7;
t_f=mu/kb;
t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
c_m = 10; c_p = 10;
a = 1; I = 0; %neural voltage model param, AVB input bias current
sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
dt=1e-3; % time step size


% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_m, c_p, t_f, t_n, t_m, a, I, sigma );

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_m, c_p, t_f, t_n, t_m, a, sigma_prime);

%display result
fig1=figure(4); clf;
subplot(3,1,1);
plot(0:dt/period:1-dt/period, Z(:,1),'-','Linewidth', 4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
ylabel('Z_\kappa'); set(gca,'FontSize',30); 
title('PRCs, Increased C_p'); 
ylim([-0.2,0.2]);
subplot(3,1,2); plot(0:dt/period:1-dt/period, Z(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:));  
set(gca,'FontSize',30); ylabel('Z_A'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.1,0.1]);
subplot(3,1,3); plot(0:dt/period:1-dt/period, Z(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:dt/period:1-dt/period, Z(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('Z_V');  xlabel('\theta'); plot([0,1],[0,0],'k:','linewidth',2);
ylim([-0.05,0.05]);
set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);


% ---- III.  CALCULATE G-FUNCTION  ----
% -- Mechanical coupling -computing weight matrix B---
gridsz=1; dim =2;
delX = 1/(gridsz*dim);
e = ones(gridsz*dim,1);
%2nd difference operator
D2 = spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%4th difference operator
D4 = D2*D2';
D4(1,1) = 7; D4(end,end) = 7; %moment-free, force-free BCs
D4 = (1/delX^4).*D4;
%mechanical coupling weight matrix B
B = inv(full(D4));

%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g_gap = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig = figure(10);
% subplot(2,1,1); 
mech3 = plot(theta,g1,'--r', 'linewidth',4); hold on;
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
plot([0,1],[0,0],'k:','linewidth',2);
% title('\tau_b = 10/3 \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25);
ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)
legend([mech1 mech2 mech3],{'tristable','increased c_m', 'increased c_p'})

%G fn - gap junctional coupling
fig2 = figure(12);
gap3 = plot(theta,g_gap/max(g_gap),'--b', 'linewidth',4); hold on; 
%find zeros & mark stability
pl_state_inds = find(abs(g_gap)/max(g_gap)<=1e-3);
stab_state_inds = find(g_gap-circshift(g_gap,-1) >= 0);
unstab_state_inds = find(g_gap-circshift(g_gap,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2);
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_g(\phi)'); set(gca,'FontSize',30)
legend([gap1 gap2 gap3],{'tristable','increased c_m', 'increased c_p'})

%G fn - proprioceptive coupling
fig3 = figure(13); load('colorblind_colormap.mat');
prop3 = plot(theta,h1p/max(h1p),'--', 'linewidth',4, 'Color', colorblind(10,:,:)); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(h1p)/max(h1p)<=1e-3);
stab_state_inds = find(h1p-circshift(h1p,-1) >= 0);
unstab_state_inds = find(h1p-circshift(h1p,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); 
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_p(\phi)'); set(gca,'FontSize',30)
legend([prop1 prop2 prop3],{'tristable','increased c_m', 'increased c_p'})

