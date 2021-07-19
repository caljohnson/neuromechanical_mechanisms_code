% ==================================================================
%
%                     full_model_no_gaps.m
%                      ------- 
%  This program runs the full model for the chain of 
%  dim  neuromechanical oscillator modules and loops over eps_p
%  First, fits t_m in [50,250] ms to match frequency and
%  then loops over eps_prop at low viscosity



% ==================================================================

addpath('./src');
clear


% -- NM MODEL PARAMETERS --
mu = 1.3e-7; %N (mm)^2 s
kb = 2.6e-7; %N (mm)^2
% mu = 10^(-9.3); %N (mm)^2 s
% kb = 9.7724e-10; %N (mm)^2
t_f=mu/kb;

t_n = 1e-2; t_m = 1e-1; %timescales for neural (10ms) and muscle activity (50-200ms)
c_ma = 5;%[5.2; 5; 4.8; 4.6; 4.4; 4;]; %musc. activity feedback strength
c_prop = 1;  %prop feedback strength (arbitrary)

%dont change these much
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
% sigma = @(x) x;
% sigma_prime = @(x) 1;
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 6; %chain of 6 units

%coupling params
little_gamma = [1, 10, logspace(2,3,4), 10^4, 10^4.4472];

dt = t_n/10; 

%loop to determine t_m so that the period is within .1 of 0.5 sec
period = 0;
true_period = 0.5;
t_m_low = 0.05;
t_m_hi = 0.25;
no_runs = 1; %number of times run this loop
%bisection on t_m to find in [50 ms ,250 ms]
while abs(period - true_period)>10^-1 && no_runs < 20
    tic
    t_m = (t_m_hi+t_m_low)/2;
    [ X_LC, period ] = single_oscillator_LC( dt, c_ma(1), c_prop, t_f, t_n, t_m, a, I, sigma );
    period
    no_runs = no_runs+1;
    if period < true_period
        t_m_low = t_m;
    else
        t_m_hi = t_m;
    end
end
%make init cond
phis = 0.7*ones(dim-1,1);
[ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
toc

periods = zeros(size(little_gamma));
wavelns = zeros(size(little_gamma));
mean_amp = zeros(size(little_gamma));
plstate = zeros(size(little_gamma,2),dim-1);

%loop over eps_prop to find wavelengths
%wavelength is within 0.1 of 1.5 wavelengths/bodylength
jj=1;
CN = (3.4*1e-9)*little_gamma(jj);
eps_props = logspace(-5, 1, 10);
eps_gap = 0;

for jj = 1:10
    tic
    eps_prop = eps_props(jj)
    %make oderhs
    oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
                t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
    %run model
    TF = 3e2;
    %run for TF time, if not enough cycles for phase differences, run again
    notlongenough=true;
    no_takes = 1; %number of TF runs
    while notlongenough
        tic;
        if no_takes==1
            [t,y1] = ode23(oderhs,[0,TF], init_cond);
        else
            %use end of last y as start
            [t2, y2] = ode23(oderhs,[TF*(no_takes-1),TF*no_takes], y1(end,:));
            %stack em together
            t = [t(1:end-1);t2;];
            y1 = [y1(1:end-1,:); y2;];
        end
        toc;
  
        tic;
        %interpolate and compute phase differences
        dt = 1e-3;
        t0 = 0:dt:TF*no_takes;
        y = interp1(t,y1,t0);
    
        try
            %compute phase differences relative to osc 1
            [phis_full,new_period]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
            toc;
            plstate(jj,:) = mean(phis_full(:,end-10:end),2)';
            periods(jj) = new_period*dt;
            freq = 1./periods(jj);
            wavelns(jj) = 1/((6/5)*sum(1-mod(plstate(jj,:),1)));
            wvln = wavelns(jj)
            notlongenough = false;
        catch
            no_takes=no_takes+1;
            fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
        end
    end
    no_runs = no_runs+1;
end

figure(1); clf;
semilogx(eps_props, wavelns,'ob-','MarkerSize',10,'MarkerFaceColor','b','LineWidth',4);
xlabel('proprioceptive coupling strength \epsilon_p'); ylabel('wavelength per bodylength (\lambda / L)');
set(gca, 'FontSize',30); xlim([10^-5, 10^2]);