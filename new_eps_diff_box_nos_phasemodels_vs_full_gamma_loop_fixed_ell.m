% ==================================================================
%
%                     new_eps_diff_box_nos_phasemodels_vs_full_gamma_loop_fixed_ell.m
%                      ------- 
%  Runs the full model and phase model for the N-box chain of NM
%  oscillators to investigate coupling
%       with fixed module length ell


% and generates the data for figures of N-oscillator wavelengths and phase-diffs
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
delX = 1/6;

dims = [2:6 10 12 20 25];
for jj=1:size(dims,2)
    dim = dims(jj); %chain of 6 units

%range of external fluid viscosities
little_gamma = logspace(0,4.4472,10);

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

% ----------------------------------
% Run phase model for different mech coupling strengths
eps_prop = 0.0134;
eps_gap = 0.005;

options = odeset('RelTol',1e-6, 'AbsTol', 1e-6); %options for odes
phis = -0.1*ones(dim-1,1);
% phis = [0.5; 0.5; 0.5; 0.5; 0.5;];
phase_model_eq2 = zeros(size(little_gamma,2),dim-1);
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
plstate2 = zeros(size(little_gamma,2),dim-1);
for kk = 1:size(little_gamma,2)
    
    CN = (3.4*1e-9)*little_gamma(kk);

    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
        t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

    %get initial condition from phase-locked states above
    phi = mod(phase_model_eq2(kk,:),1)
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

mat_file_name = strcat(num2str(dims(jj)), 'box_phase_vs_full_gamma_loop_data_neweps.mat');
save(mat_file_name);

disp(strcat('finished run with dim = ', num2str(dims(jj))))
end


