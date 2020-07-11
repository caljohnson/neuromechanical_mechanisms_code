% ==================================================================
%
%                     diff_eps_20_box_phasemodels_vs_full_gamma_loop_fixed_ell.m
%                      ------- 
%  Runs the full model and phase model for the 20-box chain of NM 
%  oscillators (w fixed module length ell) to investigate  neural coupling effects


% and generates the data for figures of 20-oscillator wavelengths and phase-diffs
% ==================================================================
addpath('./src');
clear


load('fullmodel_eps_gamma_loop.mat')

TF_phasemodel = 1e3;
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
delX = 1/6;

dim = 20;
eps_gaps = [1e-4; 1e-3; 1e-2; 1e-1;]';
eps_props = zeros(size(eps_gaps));
% eps_props = 2.7.*eps_gaps;


%range of external fluid viscosities
% mu_f = [1, 10, logspace(2,3,4), 10^4, 10^4.4472];
mu_f = logspace(0,4.4472,5);

%store wavelengths
wavelengths_phasemodel = zeros(size(eps_gaps,2),size(mu_f,2));

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
% delX = 1/(gridsz*dim);
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


for jj = 1:size(eps_gaps,2)
    %loop over eps_gap
    eps_gap = eps_gaps(jj)
%     eps_prop = eps_props(jj)

% ----------------------------------
% Run phase model for different mech coupling strengths
phase_model_eq2 = zeros(size(mu_f,2),dim-1);


%first set eps_p to get approx right wavelength
wvln = 0;
true_wvln = 1.54;
eps_prop = 2.7*eps_gap;
eps_step = .01;
no_runs = 1;
eps_prop_low = eps_prop/10;
eps_prop_hi = eps_prop*10;
tic
while abs(wvln - true_wvln)>5*10^-2 && no_runs < 20
    if no_runs>1
        if wvln>true_wvln
            eps_prop_low = eps_prop;
        else
            eps_prop_hi = eps_prop;
        end
        eps_prop = (eps_prop_hi+eps_prop_low)/2
    end

options = odeset('RelTol',1e-6, 'AbsTol', 1e-6); %options for odes
phis = -0.1*ones(dim-1,1);
CN = (3.4*1e-9)*mu_f(1);

%create phase model eqns
phasemodel = coupled_oscillator_phase_difference_odes(dim,H_m,H_p,H_n,...
        mu,CN,eps_prop/t_n,eps_gap/t_n);
    
%run phase model
[t,y] = ode23(phasemodel, [0 TF_phasemodel], phis, options);
    
%compute phase-locked state
phase_model_eq2(1,:) = y(end,:);
phis = phase_model_eq2(1,:);
wvln= 1/6./(sum(1-mod(phis,1),2)./(dim-1));
toc
no_runs =  no_runs +1;
end
%now rest of fluid viscosities
for kk = 2:size(mu_f,2)
    tic
    CN = (3.4*1e-9)*mu_f(kk);

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
% 
% %%---- Find same behavior in full model
% plstate2 = zeros(size(mu_f,2),dim-1);
% for kk = 1:size(mu_f,2)
%     
%     CN = (3.4*1e-9)*mu_f(kk);
% 
%     %make oderhs
%     oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
%         t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
% 
%     %get initial condition from phase-locked states above
%     phi = mod(phase_model_eq2(kk,:),1)
%     init_cond  = phases_to_init_cond( dim, phi, X, period, dt, gridsz );
%     
%     %run model
%     tic;
%     [t,y] = ode23(oderhs,[0,TF_full], init_cond);
%     
%     t0 = 0:1e-3:TF_full;
%     y = interp1(t,y,t0);
%     
%     %compute phis
%     phis_full  = full_timetrace_to_phasediffs( dim, y, gridsz,10 );
%     toc;
%     %set last phi as phase locked state
%     plstate2(kk,:) = mean(phis_full(:,end-10:end),2)';
% 
% end
eps_props(jj) = eps_prop;
wavelengths_phasemodel(jj,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
% wavelengths_fullmodel(jj,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
end


% inds = [2 4 6 9];
legendCell = cellstr(strcat('\',num2str(eps_gaps', 'epsilon_g=%.4f')));
legendCell{5} = 'FY 2010';
load('colorblind_colormap.mat');
colors = colorblind([1 2 6 7],:,:);
figure(1); clf;
h = semilogx(mu_f,wavelengths_phasemodel,'LineWidth',4.0,'Marker','o','Markersize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
% h2 = semilogx(mu_f,wavelengths_fullmodel,'+','LineWidth',2,'MarkerSize',10); 
% set(h2, {'color'}, num2cell(colors,2));
%plot FY data
semilogx([1, 10, logspace(2,3,4), 10^4, 10^4.4472]...
    , [1.54, 1.375, 1.25, 1.2, 1, .9, .8, .75],':k','LineWidth',3,'Marker','x','Markersize',20);
set(gca,'FontSize',30);
legend(legendCell,'location','northeastoutside','FontSize',30);
xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
ylabel({['wavelength'] ['(\lambda / L)']});


