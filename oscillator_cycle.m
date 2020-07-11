% ==================================================================
%
%                     oscillator_cycle.m
%
%  This program plots the oscillator cycle curves and PRCs
%
% ==================================================================

addpath('./src');
clear

% -- NM MODEL PARAMETERS --
t_f=.05; t_n = 1e-2; t_m = .1; %timescales for length, neural, and muscule activity
c_ma = 5; c_prop = 1;  %musc. activity feedback strength, prop feedback strength
a = 1; I = 0; %neural voltage model param, AVB input    bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
% sigma = @(x) x; %thresholding function for muscle activity into tension
% sigma_prime = @(x) 1; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 1; %chain of 6 units

% time step size
dt=1e-3;

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );

%display result
figure(1); clf;
subplot(3,1,1); plot(dt:dt:period, X(:,1),'Linewidth', 4); ylabel('\kappa');
xlim([0, period]);
set(gca,'FontSize',30)
subplot(3,1,2); plot(dt:dt:period, X(:,2), 'g','Linewidth', 4); hold on;
plot(dt:dt:period, X(:,3), 'r','Linewidth', 4); 
xlim([0, period]);
set(gca,'FontSize',30)
ylabel('A'); legend('V', 'D');
subplot(3,1,3); plot(dt:dt:period, X(:,4), 'g','Linewidth', 4); hold on;
plot(dt:dt:period, X(:,5), 'r','Linewidth', 4); 
ylabel('V'); legend('V', 'D');
xlim([0, period]);
% subplot(4,1,4); plot(sigma(X(:,2)), 'g','Linewidth', 4); hold on; plot(sigma(X(:,3)), 'r','Linewidth', 4); 
% ylabel('\sigma(A)'); legend('V', 'D');
suptitle(strcat('Oscillator Cycle, c_{MA} = ', num2str(c_ma),', T = ', num2str(period)));
set(gca,'FontSize',30)

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

%display result
figure(2); clf;
subplot(3,1,1); plot(dt:dt:period, Z(:,1), 'Linewidth', 4); ylabel('\kappa');
xlim([0, period]); set(gca,'FontSize',30);
subplot(3,1,2); plot(dt:dt:period, Z(:,2), 'g', 'Linewidth', 4); hold on; 
plot(dt:dt:period, Z(:,3), 'r','Linewidth', 4); set(gca,'FontSize',30)
ylabel('A'); legend('V', 'D'); xlim([0, period]);
subplot(3,1,3); plot(dt:dt:period, Z(:,4), 'g','Linewidth', 4); hold on; 
plot(dt:dt:period, Z(:,5), 'r','Linewidth', 4); set(gca,'FontSize',30)
ylabel('V'); legend('V', 'D'); xlim([0, period]);
% set(gcf, 'Position', get(0, 'Screensize'));
suptitle(strcat('Phase Response Curves, c_{MA} = ', num2str(c_ma), ', T = ', num2str(period)));

