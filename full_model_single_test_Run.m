% ==================================================================
%
%                     full_model_single_test_Run.m
%                      ------- 
%  This program runs the full model for the chain of 
%  dim  neuromechanical oscillator modules once. that's it

% ==================================================================

addpath('./src');
clear

% -- NM MODEL PARAMETERS --

mu = 2e-9; %N (mm)^2 s
kb = 4e-9; %N (mm)^2
t_f=mu/kb 

t_n = 1e-2; t_m = 5e-2; %timescales for neural (10ms) and muscle activity (50-200ms)
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
eps_prop = .1725;
eps_gap = .05;

%make init cond
dt = t_n/10; tic; display('finding limit cycle:');
[ X_LC, period ] = single_oscillator_LC( dt, c_ma(1), c_prop, t_f, t_n, t_m, a, I, sigma );
period
%make init cond
phis = 0.7*ones(dim-1,1);
[ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
toc

%Gamma = Cn/kb = 3.4gamma/kb*(conversionfactor)
CN = (3.4*1e-9)*little_gamma(end);
%make oderhs
oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
%run ode
TF = 1e3;
tic; display('running ode solver:');
[t,y] = ode23(oderhs,[0,TF], init_cond);
toc;
%visual curvature kymographs
[T,D] = meshgrid(1:dim+1,t);
Kappa = y(:,1:gridsz*dim);
% mean_amp(jj) = mean(max((abs(Kappa))));
Kappa(:,gridsz*dim+1) = Kappa(:,gridsz*dim);
figure(3); clf; surf(T,D,Kappa);
view(2); shading flat; colormap(blueblackred); colorbar();
ylim([t(end) - 3*period, t(end)])
% pause()


%interpolate and compute phase differences
dt = 1e-3; tic; display('finding phase differences:');
t0 = 0:dt:TF;
y2 = interp1(t,y,t0);
%compute phase differences relative to osc 1
[phis_full,period]  = full_timetrace_to_phasediffs( dim, y2, gridsz,20 );
figure(4); clf; plot(phis_full', 'o');
toc;

plstate = mean(phis_full(:,end-10:end),2)';
waveln = 1/((6/5)*sum(1-mod(plstate(:),1)))
period = period*dt
freq = 1./period