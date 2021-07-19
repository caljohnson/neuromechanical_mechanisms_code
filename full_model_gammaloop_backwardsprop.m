% ==================================================================
%
%                     full_model_gammaloop_backwardsprop.m
%                      ------- 
%  This program runs the full model with backwards proprioception for the chain of 
%  dim  neuromechanical oscillator modules and loops over external
%  viscosity gamma.  First, fits t_m in [50,250] ms to match frequency and
%  then fits eps_prop to match wavelength in water
%   Then loops over external viscosity gamma


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
phis = 0.6*ones(dim-1,1);
[ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
toc

periods = zeros(size(little_gamma));
wavelns = zeros(size(little_gamma));
mean_amp = zeros(size(little_gamma));
plstate = zeros(size(little_gamma,2),dim-1);

%first loop to determine eps_prop and eps_gap so that the low-viscosity
%wavelength is within 0.1 of 1.5 wavelengths/bodylength
jj=1;
CN = (3.4*1e-9)*little_gamma(jj);
eps_prop = 0.06;
eps_gap = 0.01;

wvln = 0;
true_wvln = 1.54;
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
    oderhs = full_state_model_odes_backwardsprop(dim, gridsz, c_ma, c_prop, mu,kb, ...
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

%show weak-coupling predictions
N_cycles = 20;
dim = 2; %chain of 2 units
% time step size
dt=1e-3;

% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X_LC, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);

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
[ H_m, H_p, H_g ] = oscillator_coupling_fns( X_LC, Z ,dt,period);
toc

%---------------------- Loop with Phase Model
theta = 0:.0001:1;
g_mech = -(H_m(-theta)-H_m(theta));
g_gap = 2*(H_g(-theta) - H_g(theta));
h1p = 2*(- H_p(theta));

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

%-- loop over fluid viscosities
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
dim = 6;
%loop over rest of gammas to show gait adaptation
for jj = 2:size(little_gamma,2)
    CN = (3.4*1e-9)*little_gamma(jj);
%make oderhs
oderhs = full_state_model_odes_backwardsprop(dim, gridsz, c_ma, c_prop, mu,kb, ...
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
%         figure(4); clf; plot(phis_full', 'o');
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

str3 = ['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
fig7 = figure(7); clf; 
subplot(2,1,1); semilogx(little_gamma, 1-state2, 'bo',...
    'Markersize',20,'LineWidth',4.0); title(str3);
hold on; semilogx(little_gamma, plstate, 'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r');
% xlabel('external fluid viscosity \mu_f (mPa s)'); 
ylabel({'phase', 'difference \phi'});  ylim([0.5, 1]);
set(gca,'FontSize',30)
subplot(2,1,2); semilogx(little_gamma,wvlns2,'bo',...
    'Markersize',20,'LineWidth',4.0); hold on;
xlabel('external fluid viscosity \mu_f (mPa s)'); 
ylabel({'wavelength','(\lambda/L)'});
semilogx(little_gamma,wavelns,'rs',...
    'Markersize',15, 'MarkerFaceColor', 'r'); ylim([0,2]);
lgd = legend('2-osc. phase model', '2-module NM model');
lgd.Position = [ 0.5644    0.4358    0.3400    0.1060];
set(gca,'FontSize',30); set(gcf,'Position', [1    59   853   646]);


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

