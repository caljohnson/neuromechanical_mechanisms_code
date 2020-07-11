% ==================================================================
%
%       is_coupling_weak_single_module_vs_fully_coupled_module.m
%                      ------- 
%  This program answers the question: Is coupling weak?
%       by comparing the single-module oscillation with a single-module
%       from the fully-coupled model

%  Plots the three separate oscillations figure: Single Isolated Oscillator, 
%   Single Oscillator (non-isolated) in NM Model at low and high viscosity
% ==================================================================

addpath('./src');
clear
load('colorblind_colormap.mat');

% -- NM MODEL PARAMETERS --

mu = 1e-7; %N (mm)^2 s
kb = 2e-7; %N (mm)^2
t_f=mu/kb 

t_n = 1e-2; t_m = 1.5e-1; %timescales for neural (10ms) and muscle activity (50-200ms)
c_ma = 5; %musc. activity feedback strength
c_prop = 1;  %prop feedback strength (arbitrary)

%dont change these much - nonlinear threshold & bistable function parameters
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 6; %chain of 6 units

%coupling params
little_gamma = [1, 10, logspace(2,3,4), 10^4, 10^4.4472]; %mPa s
eps_prop = .05;
eps_gap = .0134;

%make init cond
dt = t_n/10; tic
[ X_LC, period ] = single_oscillator_LC( dt, c_ma(1), c_prop, t_f, t_n, t_m, a, I, sigma );
X = X_LC;

%plot single-oscillator LC
fig1=figure(1); clf;
subplot(3,1,1);
plot(dt/period:dt/period:1, X(:,1),'-','Linewidth', 4); hold on;
str1 = {['Isolated Oscillation,'] ['T = ', num2str(period), 's']};
ylabel('\kappa'); title(str1);
xlim([0, 1]);
set(gca,'FontSize',30)
subplot(3,1,2); plot(dt/period:dt/period:1, X(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(dt/period:dt/period:1, X(:,3),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
xlim([0, 1]);
set(gca,'FontSize',30)
ylabel('A'); %lgd1= legend('Ventral', 'Dorsal');
% lgd1.Position = [ 0.6688    0.5588    0.1953    0.1060];
subplot(3,1,3); plot(dt/period:dt/period:1, X(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(dt/period:dt/period:1, X(:,5),'-','Linewidth', 4,'Color', colorblind(2,:,:)); 
ylabel('V'); %lgd2 = legend('Ventral', 'Dorsal');
% lgd2.Position = [0.6688     0.2585    0.1953    0.1060];
xlim([0, 1]); xlabel('\theta');
set(gca,'FontSize',30);set(gcf,'Position',[1    59   640   646]);
saveas(fig1,'figs/single_isol_osc_cycle.png');


%make init cond
phis = 0.7*ones(dim-1,1);
[ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
toc
%------------fully-coupled model, gamma lowest------------------------
CN = (3.4*1e-9)*little_gamma(1);
%make oderhs
oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop,  mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
%run ode
TF = 3e2;
tic;
[t,y] = ode23(oderhs,[0,TF], init_cond);
toc;

%interpolate and compute phase differences
dt = 1e-3; tic;
t0 = 0:dt:TF;
y = interp1(t,y,t0);
%compute phase differences relative to osc 1
[~,period_new]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
% figure(4); clf; plot(phis_full', 'o');
toc;

%get X_LC from fully-coupled module
X_new = y(end-period_new:end,[1,dim+1,2*dim+1,3*dim+1,4*dim+1]);
%shift so that it starts at the max
[~, max_ind] = max(X_new(:,1));
X_new = circshift(X_new,-max_ind,1);

%compare periods
period
period_new*dt

%plot single-oscillator LC from fully-coupled model, gamma lowest
fig2=figure(2); clf;
% figure(1);
subplot(3,1,1); 
plot(0:1/period_new:1, X_new(:,1),'-','Linewidth',4); 
str2 = {['Oscillation from NM Model,'] ['T = '...
    num2str(period_new*dt) 's , low viscosity']};% \mu_f = ' ...
%     num2str(little_gamma(1)) 'mPa s']};
ylabel('\kappa'); title(str2);
xlim([0, 1]);
set(gca,'FontSize',30)
subplot(3,1,2); plot(0:1/period_new:1, X_new(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:1/period_new:1, X_new(:,3), '-','Linewidth', 4, 'Color', colorblind(2,:,:)); 
xlim([0, 1]);
set(gca,'FontSize',30)
ylabel('A');%lgd1= legend('Ventral', 'Dorsal');
% lgd1.Position = [ 0.6688    0.5588    0.1953    0.1060];
subplot(3,1,3); plot(0:1/period_new:1, X_new(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:1/period_new:1, X_new(:,5), '-','Linewidth', 4, 'Color', colorblind(2,:,:)); 
ylabel('V');%lgd2 = legend('Ventral', 'Dorsal');
% lgd2.Position = [0.6688     0.2585    0.1953    0.1060];
xlim([0, 1]); xlabel('\theta');
set(gca,'FontSize',30);set(gcf,'Position',[1    59   640   646]);
saveas(fig2,'figs/single_osc_cycle_fullmodel_muf1.png');

%------------fully-coupled model, gamma highest------------------------
%Gamma = Cn/kb = 3.4gamma/kb*(conversionfactor)
CN = (3.4*1e-9)*little_gamma(4);
%make oderhs
oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop,  mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
%run ode
TF = 1e3;
tic;
[t,y] = ode23(oderhs,[0,TF], init_cond);
toc;

%interpolate and compute phase differences
dt = 1e-3; tic;
t0 = 0:dt:TF;
y = interp1(t,y,t0);
%compute phase differences relative to osc 1
[~,period_new2]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
% figure(4); clf; plot(phis_full', 'o');
toc;

%get X_LC from fully-coupled module
X_new2 = y(end-period_new2:end,[1,dim+1,2*dim+1,3*dim+1,4*dim+1]);
%shift so that it starts at the max
[~, max_ind] = max(X_new2(:,1));
X_new2 = circshift(X_new2,-max_ind,1);

%compare periods
period
period_new*dt
period_new2*dt

%plot single-oscillator LC from fully-coupled model, gamma highest
fig3= figure(3); clf;
subplot(3,1,1); 
plot(0:1/period_new2:1, X_new2(:,1),'-','Linewidth', 4); ylabel('\kappa');
xlim([0, 1]); 
str3 = {['Oscillation from NM Model,'] ['T = '...
    num2str(period_new2*dt) 's , high viscosity']};% \mu_f = ' ...
%     num2str(little_gamma(end),'%.2g') 'mPa s']};
title(str3);
% lgd.Location = 'northoutside'
set(gca,'FontSize',30)
subplot(3,1,2); plot(0:1/period_new2:1, X_new2(:,2), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:1/period_new2:1, X_new2(:,3), '-','Linewidth', 4, 'Color', colorblind(2,:,:)); 
xlim([0, 1]);
set(gca,'FontSize',30)
ylabel('A');%lgd1= legend('Ventral', 'Dorsal');
% lgd1.Position = [ 0.6688    0.5588    0.1953    0.1060];

subplot(3,1,3); plot(0:1/period_new2:1, X_new2(:,4), '--','Linewidth', 4,'Color', colorblind(10,:,:)); hold on;
plot(0:1/period_new2:1, X_new2(:,5), '-','Linewidth', 4, 'Color', colorblind(2,:,:)); 
ylabel('V'); %lgd2 = legend('Ventral', 'Dorsal');
% lgd2.Position = [0.6688     0.2585    0.1953    0.1060];
xlim([0, 1]); xlabel('\theta'); set(gca,'FontSize',30);
set(gcf,'Position',[1    59   640   646]);
saveas(fig3,'figs/single_osc_cycle_fullmodel_muf2e4.png');
% saveas(fig1,'figs/compare_osc_cycles_kappa');



