% ==================================================================
%
%                     twobox_invest_loop.m
%                      ------- 
%  Runs the full model and phase model for the 2-box pair of NM
%  oscillators to investigate coupling dynamics
%   first runs loop at gamma=1mPa s to determine eps_prop, eps_gap to get
%   correct wavelength in water
%   then loops over rest of viscosities to get coordination trend
%   also compares to a direct simulation of the fully-coupled model with
%   only 2 modules

%  this makes several figures:
%   G fns of each coupling modality
%   eps_prop/eps_gap vs. phase difference and vs. wavelength
%   gamma vs phase difference and vs. wavelength

% ==================================================================
addpath('./src');
clear

TF = 2e2; %max runtime for full 2-box sims
N_cycles = 10; %Number of cycles to average phase over
tic
% -- NM MODEL PARAMETERS --
mu = 1.3e-7;
kb = 2.6e-7;
t_f=mu/kb;
t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
c_ma = 5; c_prop = 1;  %musc. activity feedback strength, prop feedback strength
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 2; %chain of 2 units

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

%---------------------- Loop with Phase Model
little_gamma = logspace(0,4.4472,10);
% Gammas = [1e-3; 1e-2; 1e-1; 1; 1e1; 1e2;];
% eps_props = [0; 0.001; 0.01; 0.1; 0.5;];
% eps_gaps = [0; 0.001; 0.01; 0.1; 0.5;];

theta = 0:.0001:1;
g_mech = -(H_m(-theta)-H_m(theta));
g_gap = 2*(H_g(-theta) - H_g(theta));
h1p = -2*H_p(-theta);

fig1 = figure(1); clf;
plot(theta,g_mech/max(g_mech),'r', 'linewidth',4); hold on; 
%find zeros & mark stability
pl_state_inds = find(abs(g_mech)/max(g_mech)<=1e-3);
stab_state_inds = find(g_mech-circshift(g_mech,-1) > 0);
unstab_state_inds = find(g_mech-circshift(g_mech,-1) <= 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
plot([0,1],[0,0],'k:','linewidth',2); 
xlabel('phase \phi');  ylim([-1 1]); ylabel('G_m(\phi)'); set(gca,'FontSize',30)
saveas(fig1,'figs/gfn_mech.png');

fig2 = figure(2); clf;
plot(theta,g_gap/max(g_gap),'b', 'linewidth',4); hold on; 
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
saveas(fig2,'figs/gfn_gap.png');

fig3 = figure(3); clf; load('colorblind_colormap.mat');
plot(theta,h1p/max(h1p),'-', 'linewidth',4, 'Color', colorblind(10,:,:)); hold on;
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
saveas(fig3,'figs/gfn_prop.png');
% title(strcat('G-functions (normalized) for NM paired-oscillator model'));% c_{MA} = ', num2str(c)));

%--- loop over eps_gap 
ratio = logspace(-.5,1,15);
eps_prop = 0.05;
% 
CN = (3.4*1e-9)*little_gamma(1);
for kk = 1:size(ratio,2)
    
    eps_gap = eps_prop/ratio(kk);
    m_strength = CN/mu*B(1,2);
    g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gap*g_gap/t_n;
%     figure(2); clf;
%     plot(theta,g(:),'m-','linewidth',2); hold on;
%     plot(theta,eps_gap*g_gap(:),'b-', 'linewidth',2);
%     plot(theta,eps_prop*h1p(:),'g-','linewidth',2);
%     plot([0,1],[0,0],'k:','linewidth',2);

%     plot(theta(abs(g)/max(g)<=5e-3),0*find(abs(g)/max(g)<=5e-3),'ko','Markersize',10);
    zeroes{kk} = theta(abs(g)/max(g)<=1e-2);
%     plot(theta(abs(g_mech)/max(g_mech)<=1e-3),0*find(abs(g_mech)/max(g_mech)<=1e-3),'ko','Markersize',10);
%     plot(theta(abs(h1p)/max(h1p)<=2e-3),0*find(abs(h1p)/max(h1p)<=2e-3),'ko','Markersize',10);
%     legend('both','gap-junction','proprioceptive','"zero"');
%     xlabel('\phi');  ylabel('d \phi / d t');
%     str1=['G-function for NM paired-oscillator model; \epsilon_{prop} = ' num2str(eps_prop) ', \epsilon_{gap} = ' num2str(eps_gap) ', \Gamma = 0'];
%     title(str1);
%     set(gca,'FontSize',30)
% pause(0.5);
% pause;
end;

%turn stable phase-locked state into a "wavelength" in the 6-box case
wvlns = zeros(size(ratio,2),1);
for kk=1:size(ratio,2)
    states = uniquetol(zeroes{1,kk},1e-1);
    %second state will be the stable one, but will be a phase-advance > 0.5
    state(kk) = states(2); %turn it into a phase-delay of magnitude <0.5
    %turn state into wavelength (in 6-box case, assuming this phase diff
    %becomes a constant phase wave)
    wvlns(kk) = 1./((1-state(kk))*6);
end

% figure(3); clf; semilogx(ratio,wvlns,'o'); hold on;
% xlabel('proprioceptive to gap-junction coupling strength ratio'); ylabel('wavelength per bodylength');
% str2=['Wavelengths for 6-Box chain from NM paired-oscillator model vs. \epsilon_{prop}, with \epsilon_{gap} = ' num2str(eps_gap) ' and \Gamma = 0'];
% title(str2);
toc

%%---- Find same behavior in full model
eps_gaps = eps_prop./ratio;

for kk = 1:size(ratio,2)
    
    eps_gap = eps_prop/ratio(kk);
    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
        t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

    %get initial condition from phase-locked state above
    state(kk)
    init_cond  = phases_to_init_cond( dim, state(kk), X, period, dt, gridsz );
    
    %run model
    tic;
    [t,y] = ode23(oderhs,[0,TF], init_cond);
    
    t0 = 0:1e-3:TF;
    y = interp1(t,y,t0);
    
    %compute phis
    phis_full  = full_timetrace_to_phasediffs( dim, y, gridsz,N_cycles );
   
    toc;
    %set mean of last few phis as phase locked state
    plstate(kk) = mean(phis_full(end-10:end));

end

%turn stable phase-locked state into a "wavelength" in the 6-box case
plstate = mod(plstate,1);
wvlns_full = 1./((1-plstate)*6);

% figure(3); semilogx(ratio,wvlns_full,'o');
% legend('2-osc. phase model', '2-module NM model');
% xlabel('log_{10}\epsilon_{prop}'); ylabel('\lambda / L');
% str2=['Wavelengths for 6-Box chain from NM paired-oscillator model vs. \epsilon_{prop}, with \epsilon_{gap} = ' num2str(eps_gap) ' and \Gamma = 0'];
% title(str2);

% figure(10); clf;
% semilogx(ratio, state, 'o')
% hold on; semilogx(ratio, plstate, 'o')
% xlabel('proprioceptive to gap-junction coupling strength ratio');
% ylabel('phase difference \phi');
% legend('2-osc. phase model', '2-module NM model');

str3 = ['\epsilon_{p} = ' num2str(eps_prop) ...
   ' and \mu_f = 1 mPa\cdot s'];
fig6 = figure(6); clf;
subplot(2,1,1); semilogx(eps_gaps, state, 'bo',...
    'Markersize',20,'LineWidth',4.0); title(str3);
hold on; semilogx(eps_gaps, plstate, 'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r');
% xlabel('gap-junction coupling strength \epsilon_{g}');
ylabel({'phase', 'difference \phi'}); ylim([0.5, 1]);
% legend('2-osc. phase model', '2-module NM model'); 
xlim([eps_gaps(end), eps_gaps(1)]);
set(gca,'FontSize',30); 

subplot(2,1,2); loglog(eps_gaps,wvlns,'bo',...
    'Markersize',20,'LineWidth',4.0); hold on;
xlabel('gap-junction coupling strength \epsilon_g');
ylabel({'wavelength', '(\lambda/L)'});
loglog(eps_gaps,wvlns_full,'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r');
xlim([eps_gaps(end), eps_gaps(1)]);
ylim([10^{-1},10^1]);
lgd = legend('2-osc. phase model', '2-module NM model');
lgd.Position = [ 0.5644    0.4358    0.3400    0.1060];
set(gca,'FontSize',30); set(gcf, 'Position', [1    59   853   646]);
% saveas(fig6,'figs/2osc_phases_neural_longwavelength.png');

return
%make plot of G-fns for various eps_g
fig9=figure(9); clf;  
g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gaps(end)*g_gap/t_n;
subplot(3,1,1); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
ylabel('G(\phi)');
str1=['\epsilon_g = ' num2str(eps_gaps(end))]; title(str1);

g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gaps(10)*g_gap/t_n;
subplot(3,1,2); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
ylabel('G(\phi)');
str2=['\epsilon_g = ' num2str(eps_gaps(10))]; title(str2);

g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gaps(1)*g_gap/t_n;
subplot(3,1,3); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
xlabel('phase difference \phi');  ylabel('G(\phi)');
str3=['\epsilon_g = ' num2str(eps_gaps(1))]; title(str3);
set(gcf,'Position',[1    59   640   646]);
saveas(fig9,'figs/2osc_fullgfns_longwavelength.png');


%---------------------------------------------------------------
%do again, now looping over mech coupling strengths
eps_prop = 0.05;
eps_gap = 0.0134;

for kk = 1:size(little_gamma,2)
    CN = (3.4*1e-9)*little_gamma(kk);
    
    m_strength = CN/mu*B(1,2);
    g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gap*g_gap/t_n;
  
    zeroes2{kk} = theta(abs(g)/max(g)<=5e-3);
end;

%turn stable phase-locked state into a "wavelength" in the 6-box case
wvlns2 = zeros(size(little_gamma,2),1);
for kk=1:size(little_gamma,2)
    states = uniquetol(zeroes2{1,kk},1e-2);
    %second state will be the stable one, but will be a phase-advance > 0.5
    if size(states,2)>1
        state2(kk) = 1-states(2); %turn it into a phase-delay of magnitude <0.5
    else
        state2(kk) = 1-states(1); %turn it into a phase-delay of magnitude <0.5
    end
    %turn state into wavelength (in 6-box case, assuming this phase diff
    %becomes a constant phase wave)
    wvlns2(kk) = 1./(state2(kk)*6);
end

%and compare with full model
for kk = 1:size(little_gamma,2)
       
    CN = (3.4*1e-9)*little_gamma(kk);
    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
        t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

    %get initial condition from phase-locked state above
    init_cond  = phases_to_init_cond( dim, state2(kk), X, period, dt, gridsz );
    
    %run model
    tic;
    [t,y] = ode23(oderhs,[0,TF], init_cond);
    toc;
    t0 = 0:1e-3:TF;
    y = interp1(t,y,t0);
    
    %compute phis
    phis_full  = full_timetrace_to_phasediffs( dim, y, gridsz,N_cycles  );
    plstate2(kk) = mean(phis_full(end-10:end));
end

wvlns_full2 = 1./((1-plstate2)*6);
str3 = ['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
fig7 = figure(7); clf; 
subplot(2,1,1); semilogx(little_gamma, 1-state2, 'bo',...
    'Markersize',20,'LineWidth',4.0); title(str3);
hold on; semilogx(little_gamma, plstate2, 'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r');
% xlabel('external fluid viscosity \mu_f (mPa s)'); 
ylabel({'phase', 'difference \phi'});  ylim([0.5, 1]);
set(gca,'FontSize',30)
subplot(2,1,2); semilogx(little_gamma,wvlns2,'bo',...
    'Markersize',20,'LineWidth',4.0); hold on;
xlabel('external fluid viscosity \mu_f (mPa s)'); 
ylabel({'wavelength','(\lambda/L)'});
semilogx(little_gamma,wvlns_full2,'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r'); ylim([0,2]);
lgd = legend('2-osc. phase model', '2-module NM model');
lgd.Position = [ 0.5644    0.4358    0.3400    0.1060];
set(gca,'FontSize',30); set(gcf,'Position', [1    59   853   646]);
saveas(fig7,'figs/2osc_phases_gait_adapt.png');


fig8=figure(8); clf;  
m_strength = (3.4*1e-9)*little_gamma(1)/mu*B(1,2);
g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gap*g_gap/t_n;
subplot(3,1,1); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
ylabel('G(\phi)');
str1=['\mu_f = ' num2str(little_gamma(1))]; title(str1);

m_strength = (3.4*1e-9)*little_gamma(end-4)/mu*B(1,2);
g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gap*g_gap/t_n;
subplot(3,1,2); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
ylabel('G(\phi)');
str2=['\mu_f = 2.8 \times 10^3 mPa s']; title(str2);

m_strength = (3.4*1e-9)*little_gamma(end)/mu*B(1,2);
g=eps_prop*h1p/t_n + m_strength*g_mech + eps_gap*g_gap/t_n;
subplot(3,1,3); plot(theta,g(:),'m-','linewidth',4); hold on;
plot([0,1],[0,0],'k:','linewidth',2);
%find zeros & mark stability
pl_state_inds = find(abs(g)/max(g)<=1e-3);
stab_state_inds = find(g-circshift(g,-1) >= 0);
unstab_state_inds = find(g-circshift(g,-1) < 0);
plot(theta(intersect(pl_state_inds, unstab_state_inds)),...
    0*theta(intersect(pl_state_inds, unstab_state_inds)),'ko',...
    'Markersize',15);
plot(theta(intersect(pl_state_inds, stab_state_inds)),...
    0*theta(intersect(pl_state_inds, stab_state_inds)),'ko',...
    'Markersize',20, 'MarkerFaceColor', 'k');
zeroes2{kk} = theta(abs(g)/max(g)<=5e-3); set(gca,'FontSize',30)
xlabel('phase difference \phi');  ylabel('G(\phi)');
str3=['\mu_f = 2.8 \times 10^4 mPa s']; title(str3);
set(gcf,'Position',[1    59   640   646]);
saveas(fig8,'figs/2osc_fullgfns_gaitadapt.png');

