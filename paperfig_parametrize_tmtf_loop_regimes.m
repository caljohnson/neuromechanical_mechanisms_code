% ==================================================================
%                paperfig_parametrize_tmtf_loop_regimes.m
%
%  Generates data used in make_mukb_trend_colormap_figs.m to make figs
%                      ---------------------------------- 
%  First, loop over kb and t_f.  Then find t_m in 50-200ms to get correct freq in water
%   then set coupling params to get correct wavelength in water
%   then loop over external viscosities and examine coordination
%       
% ==================================================================

addpath('./src');
clear
% load('paper_trials_data.mat');
trial_number = 0;%trial(end).number +1;
% trial(trial_number).number = trial_number;

% -- NM MODEL PARAMETERS --
mus = logspace(-9.3010, log10(1.3e-7), 6)';
% kbs = logspace(log10(7.53e-10), log10(2.6e-6), 6)';
t_fs = logspace(-2,1,15)';

for mm = 1:size(t_fs,1)
    if trial_number ~= 0 
        trial_number = trial(end).number +1;
    else
        trial_number = 1
    end 
    trial(trial_number).number = trial_number

kbs = mus./t_fs(mm); %fixed so that t_f is fixed always

t_n = 1e-2; %timescales for neural (10ms)
t_m = 1e-1; %DEFAULT timescale for muscle activity (50-200ms)
c_ma = 5;%musc. activity feedback strength
c_prop = 1;  %prop feedback strength (arbitrary)

%dont change these much - nonlinear fn and threshold params
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 6; %chain of 6 units

%external fluid viscosities
mu_f = [1, 10, logspace(2,3,4), 10^4, 10^4.4472];

%make init cond
dt = t_n/10; 

%loop over mus
for kk = 1:size(mus,1)
        kb = kbs(kk)
        mu = mus(kk)
        t_f=mu/kb
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

periods = zeros(size(mu_f));
wavelns = zeros(size(mu_f));
mean_amp = zeros(size(mu_f));
plstate = zeros(size(mu_f,2),dim-1);

%first loop to determine eps_prop and eps_gap so that the low-viscosity
%wavelength is within 0.1 of 1.5 wavelengths/bodylength
jj=1;
CN = (3.4*1e-9)*mu_f(jj);

wvln = 0;
true_wvln = 1.54;
eps_prop = 0.05;
eps_gap = 0.0134;
eps_step = .01;
no_runs = 1;
while abs(wvln - true_wvln)>5*10^-2 && no_runs < 20
    if no_runs>1
        if wvln>true_wvln
            eps_prop = eps_prop+eps_step
            eps_step = eps_step*.9
        else
            eps_prop = abs(eps_prop-eps_step)
            eps_step = eps_step*.9
        end
    end

    tic
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

%loop over rest of gammas to show gait adaptation
for jj = 2:size(mu_f,2)
    CN = (3.4*1e-9)*mu_f(jj);
%make oderhs
oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

%run model
TF = 1e2;
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
        notlongenough = false;
    catch
        no_takes=no_takes+1;
        if no_takes >= 10
            fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Stopping simulation. ']);
            break
        else
            fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
        end
    end
end

end

%----
%decide whether coordination is coherent or not
%check wavelength trend
no_trend = 0;
eps_wvln = 0.01;
for jj=1:size(wavelns,2)-1
    if wavelns(jj)<wavelns(jj+1)-eps_wvln;
        no_trend = 1;
        break
    end
end
%also check phase differences for spread - indicating no traveling wave
for jj=1:size(wavelns,2)-1
    if max(plstate(jj,:)) - min(plstate(jj,:)) >= 0.4
        no_trend = 1;
        break
    end
end

%   decide whether or not 2-box phase model agrees
[phase_agrees,phase_state] = ...
            does_2box_agree(mu,kb, t_m, no_trend);
         
%save information about this run
trial(trial_number).run(kk).number = kk;
trial(trial_number).run(kk).mu = mu;
trial(trial_number).run(kk).k_b = kb;
trial(trial_number).run(kk).tau_f = t_f;
trial(trial_number).run(kk).tau_m = t_m;
trial(trial_number).run(kk).single_osc_period = period;
trial(trial_number).run(kk).eps_prop = eps_prop;
trial(trial_number).run(kk).eps_gap = eps_gap;
trial(trial_number).run(kk).frequencies = 1./periods;
trial(trial_number).run(kk).wavelengths = wavelns;
trial(trial_number).run(kk).ext_gamma = mu_f;
trial(trial_number).run(kk).phase_diffs = 1-mod(plstate,1);
%Conclusion: no trend = 1 for no wavelength trend
trial(trial_number).run(kk).no_trend = no_trend; 
trial(trial_number).run(kk).phase_agrees = phase_agrees;
trial(trial_number).run(kk).phase_state = phase_state; 

end
end

fprintf('finished all trials');
save('paper_trials_data_4.mat', 'trial')