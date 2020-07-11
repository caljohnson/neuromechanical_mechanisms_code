% ==================================================================
%
%                    compare_gfunctions.m
%                      ------- 
%  Compares the g-functions for a 2-osc pair
%  in different timescale orderings

%  Makes several figures:
%       Full G-fn in three parameter regimes
%       Bif. diag of phase-locked states vs. tau_b (aka tau_f)
% ==================================================================
addpath('./src');
clear

% ------ tau_f > tau_m ----------------
tic
% -- NM MODEL PARAMETERS --
mu = 1.3e-7; %my larger mu so that mu/kb = 0.5
kb = 2.6e-7;
t_f=mu/kb 

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
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

% ---- III.  CALCULATE G-FUNCTION  ----
% -- Mechanical coupling -computing weight matrix B---
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
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig = figure(1); clf;
% subplot(2,1,1); 
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
plot([0,1],[0,0],'k:','linewidth',2);
title('\tau_b = 10/3 \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25);
ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)
saveas(fig,'figs/compare_gfns_tf_tm_1.png');
 
% ------ tau_f < tau_m ----------------
% ----------------------------------
t_f = 0.05; %t_f < t_m=0.15

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig2 = figure(2); clf;
% subplot(2,1,2); 
plot(theta,g1,'r', 'linewidth',4); hold on;
%find zeros & mark stability
pl_state_inds = find(abs(g1)/max(g1)<=1e-3);
stab_state_inds = find(g1-circshift(g1,-1) >= 0);
unstab_state_inds = find(g1-circshift(g1,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); hold off;
title('\tau_b = 1/3 \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)
saveas(fig2,'figs/compare_gfns_tf_tm_2.png');

% ------ tau_f \approx tau_m ----------------
% ----------------------------------
t_f = 0.14; %t_f \approx t_m=0.15

% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
toc

theta = 0:.0001:1;
g1 = -(H_m(-theta)-H_m(theta));
g2 = 2/t_n*(H_g(-theta) - H_g(theta));
h1p = -2/t_n*H_p(-theta);

fig3 = figure(3); clf;
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
title('\tau_b \approx \tau_m', 'FontSize',25);
xlabel('phase difference \phi', 'FontSize', 25); ylabel('G_{m}(\phi)','FontSize', 25);
set(gca,'FontSize',25)
saveas(fig3,'figs/compare_gfns_tf_tm_3.png');


% ------ tau_f vs tau_m bif diag ----------------
t_fs = linspace(0.05,0.5,100);

for jj=1:size(t_fs,2)
    t_f = t_fs(jj);
    fprintf(['t_f = ', num2str(t_f)]);
    tic
    % ----  I. FIND PERIODIC ORBIT  ----
    [ X, period ] = single_oscillator_LC( ...
        dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );
    toc
    fprintf('found periodic orbit');

    % ----  II.  CALCULATE iPRC ---- 
    tic
    Z = single_oscillator_PRC(X, dt,...
        c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);
    toc
    fprintf('found iPRC');
    % ---- III.  CALCULATE G-FUNCTION  ----
    %compute coupling functions
    [ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period);
    theta = 0:.0001:1;
    g1 = -(H_m(-theta)-H_m(theta));
    g2 = 2/t_n*(H_g(-theta) - H_g(theta));
    h1p = -2/t_n*H_p(-theta);

    zero_inds = find(abs(g1)/max(g1)<=1e-3);
    g1 = [g1(1:end-1);g1]; %double length to make indexing easier
    %find stable states from here
    bif_diag_states(jj).stable_states = [];
    bif_diag_states(jj).unstable_states = [];
    for kk = 1:size(zero_inds)
        %check sign of G(phi) through this zero
        if g1(zero_inds(kk))>g1(zero_inds(kk)+1)
            bif_diag_states(jj).stable_states = ...
                [bif_diag_states(jj).stable_states, theta(zero_inds(kk))];
        else
            bif_diag_states(jj).unstable_states = ...
                [bif_diag_states(jj).unstable_states, theta(zero_inds(kk))];
        end
    end
end
fig4 = figure(4); clf; hold on;
for jj = 1:size(t_fs,2);
    N_s = size(bif_diag_states(jj).stable_states,2);
    N_us = size(bif_diag_states(jj).unstable_states,2);
    scatter(t_fs(jj).*ones(1,N_s),bif_diag_states(jj).stable_states,...
        100,'bo', 'filled');
     scatter(t_fs(jj).*ones(1,N_us),bif_diag_states(jj).unstable_states,...
        100,'ro');
end
plot(t_m.*ones(10,1),linspace(0,1,10),'k--', 'LineWidth',4.0);
xlabel('\tau_b (sec)'); ylabel('phase locked state \phi^*'); set(gca,'FontSize',25);
xlim([0.05,0.5]);
tick_locs = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5];
tick_labs = {0.1, '\tau_m', 0.2, 0.3, 0.4, 0.5};
set(gca, 'XTick', tick_locs, 'xticklabel', tick_labs);
legend('stable','unstable');
saveas(fig4,'figs/bif_diag_mech_coupling_phi_vs_tf_tm150ms.png');
