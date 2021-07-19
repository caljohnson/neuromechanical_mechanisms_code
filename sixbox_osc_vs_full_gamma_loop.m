% ==================================================================
%
%                     sixbox_osc_vs_full_gamma_loop.m
%                      ------- 
%  Runs the full model and phase model for the 6-box chain of NM
%  oscillators to investigate coupling


% and generates the figure of 6-oscillator wavelengths and phase-diffs
% ==================================================================
addpath('./src');
clear

TF_phasemodel = 2e3;
TF_full = 5e2;

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
dim = 6; %chain of 6 units

%range of external fluid viscosities
little_gamma = logspace(0,4.4472,15);

% time step size
dt=1e-3;
tic
% ----  I. FIND PERIODIC ORBIT  ----
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );


% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

% ---- III.  CALCULATE G-FUNCTION  ----

% -- Mechanical coupling -computing weight matrix B---
%second-difference matrix A
delX = 1/(gridsz*dim);
e = ones(gridsz*dim,1);
% A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], dim, dim+2);
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
ey = speye(gridsz*dim,gridsz*dim);
%spread operator
rowinds = 1:gridsz*dim;
colinds = repelem(1:dim, gridsz);
npdf_ids = -1:2/(gridsz-1):1;
S = sparse(rowinds, colinds, repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));
%proprioceptive averaging matrix
rowinds = repelem(1:dim,gridsz);
colinds = 1:gridsz*dim;
% Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));
Pa = sparse(rowinds, colinds, (1/gridsz).*repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));
%mechanical coupling weight matrix B
% B = Pa*((A*A')\S);
B = inv(full(A*A'));

%compute coupling functions
[ H_m, H_p, H_n ] = oscillator_coupling_fns( X, Z ,dt,period);
toc
% ----------------------------------
% Run phase model for different epsgap vs eps_prop
ratio = logspace(-0.2553,1,10);
eps_prop = 0.05;
% eps_gaps = linspace(0.05, 0.9, 15);
CN = (3.4*1e-9)*little_gamma(1);
% Gamma = (3.4/(mu*1e9))*little_gamma(1);
options = odeset('RelTol',1e-6, 'AbsTol', 1e-6); %options for odes
phis = -[0.1; 0.1; 0.1; 0.1; 0.1;];
phase_model_eq = zeros(size(ratio,2),5);
for kk = 1:size(ratio,2)
    tic
    eps_gap = eps_prop/ratio(kk);

    %create phase model eqns
    phasemodel = coupled_oscillator_phase_difference_odes(dim,H_m,H_p,H_n,...
        mu,CN,eps_prop/t_n,eps_gap/t_n);
    
    %run phase model
    [t,y] = ode23(phasemodel, [0 TF_phasemodel], phis, options);
    
%     figure(1); clf; plot(y); pause(1);
    %compute phase-locked state
    phase_model_eq(kk,:) = y(end,:);
    phis = phase_model_eq(kk,:);
    toc
end;
figure(1); clf;
subplot(2,1,1); semilogx(ratio,1-mod(phase_model_eq,1),'o'); hold on;
xlabel('proprioceptive to gap-junction coupling strength ratio'); ylabel('pair-wise phase difference');
% str1=['Phase Differences for 6-Box Chain vs. \epsilon_{prop}, with \epsilon_{gap} = ' num2str(eps_gap) ' and \Gamma = 0'];
% str1 = ['\epsilon_{prop}, with \epsilon_{gap} = ' num2str(eps_gap) ' and \Gamma = 0'];
% title(str1);

%turn into wavelength
wvln = 1./((6/5)*sum(1-mod(phase_model_eq,1),2));
figure(1); subplot(2,1,2); semilogx(ratio,wvln,'o'); hold on;
xlabel('proprioceptive to gap-junction coupling strength ratio'); ylabel('wavelength');
% str2=['Wavelengths for 6-Box chain vs. \epsilon_{prop}, with \epsilon_{gap} = ' num2str(eps_gap) ' and \Gamma = 0'];
% title(str2);
% semilogx(eps_prop./eps_gaps, 1.4661)

%%---- Find same behavior in full model
plstate = zeros(size(ratio,2),5);
for kk = 1:size(ratio,2)
    
    eps_gap = eps_prop/ratio(kk);
    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
        t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

    %get initial condition from phase-locked states above
    phi = 1-mod(phase_model_eq(kk,:),1)
    init_cond  = phases_to_init_cond( dim, phi, X, period, dt, gridsz );
    
    %run model
    tic;
    [t,y] = ode23(oderhs,[0,TF_full], init_cond);
    
    t0 = 0:1e-3:TF_full;
    y = interp1(t,y,t0);
    
    %compute phis
    phis_full  = full_timetrace_to_phasediffs( dim, y, gridsz,10);
%     figure(2); clf; plot(phis_full', 'o'); pause(1);
    toc;
    %set last phi as phase locked state
    plstate(kk,:) = mean(phis_full(:,end-10:end),2)'

end

%add to phase plot
figure(1); subplot(2,1,1);
semilogx(ratio,1-mod(plstate,1),'o');

%turn into wavelength & add to plot
wvlns_full = 1./sum(1-mod(plstate,1),2);
figure(1); subplot(2,1,2); semilogx(ratio, wvlns_full,'o');


% ----------------------------------
% Run phase model for different mech coupling strengths
eps_prop = 0.05;
eps_gap = 0.017;

options = odeset('RelTol',1e-6, 'AbsTol', 1e-6); %options for odes
phis = -[0.1; 0.1; 0.1; 0.1; 0.1;];
% phis = [0.5; 0.5; 0.5; 0.5; 0.5;];
phase_model_eq2 = zeros(size(little_gamma,2),5);
for kk = 1:size(little_gamma,2)
    tic
    CN = (3.4*1e-9)*little_gamma(kk);

    %create phase model eqns
    phasemodel = coupled_oscillator_phase_difference_odes(dim,H_m,H_p,H_n,...
        mu,CN,eps_prop/t_n,eps_gap/t_n);
    
    %run phase model
    [t,y] = ode23(phasemodel, [0 TF_phasemodel], phis, options);
    
    %compute phase-locked state
    phase_model_eq2(kk,:) = y(end,:);
    phis = phase_model_eq2(kk,:);
    toc
end;

%%---- Find same behavior in full model
plstate2 = zeros(size(little_gamma,2),5);
for kk = 1:size(little_gamma,2)
    
    CN = (3.4*1e-9)*little_gamma(kk);

    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
        t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

    %get initial condition from phase-locked states above
    phi = 1-mod(phase_model_eq2(kk,:),1)
    init_cond  = phases_to_init_cond( dim, phi, X, period, dt, gridsz );
    
    %run model
    tic;
    [t,y] = ode23(oderhs,[0,TF_full], init_cond);
    
    t0 = 0:1e-3:TF_full;
    y = interp1(t,y,t0);
    
    %compute phis
    phis_full  = full_timetrace_to_phasediffs( dim, y, gridsz,10 );
    toc;
    %set last phi as phase locked state
    plstate2(kk,:) = mean(phis_full(:,end-10:end),2)';

end

save('sixbox_osc_vs_full_gamma_loop_data.mat');

load('sixbox_osc_vs_full_gamma_loop_data.mat');
load('colorblind_colormap.mat')
fig=figure(2); clf;
colors = colorblind([1 2 6 7 8],:,:);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',4,'MarkerSize',20); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'s','LineWidth',2,'MarkerSize',15);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f (mPa s)'); 
set(gca,'FontSize',30);
text(3.5*10^4,0.725,'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
text(3.5*10^4,0.675,'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
text(3.5*10^4,0.775,'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
text(3.5*10^4,0.625,'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));


%turn into wavelength
wvln2 = 1./((6/5)*sum(1-mod(phase_model_eq2,1),2));
figure(1); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
semilogx(little_gamma,wvln2,'o-','LineWidth',4,'MarkerSize',20); hold on;
wvlns_full2 = 1./((6/5)*sum(1-mod(plstate2,1),2));
semilogx(little_gamma,wvlns_full2,'s','LineWidth',4,'MarkerSize',15);
xlabel('External Fluid Viscosity \mu_f (mPa s)'); 
ylabel({['wavelengths'], ['(\lambda/L)']});
set(gca,'FontSize',30); %title(str1);
% lgd = legend('Phase Model', 'Neuromechanical Model');
% lgd.Location = 'northoutside';
% lgd.FontSize = 25;
% saveas(fig,'figs/6oscphase_vs_fullNM_wvln_vs_gamma.png');
