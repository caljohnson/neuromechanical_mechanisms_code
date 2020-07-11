function [ yes_it_agrees,stable_state ] = does_2box_agree( mu,kb,t_m, no_trend )
%does_2box_agree Check whether 2-box phase reduction gives the same coordination
%prediction as the full simulation
%Inputs: mu, kb, t_m - oscillator parameters
%        no_trend = 1 if trend is "incorrect", i.e., if increasing
%        mechanical coupling strength does NOT decrease wavelength
%Outputs: yes_it_agrees = 1 if mech. coupling gives stable antiphase coord
%         in the 2box case AND coord_trend = 0
%         yes_it_agrees = 0 if mech. coupling gives stable antiphase coord
%         in the 2box case AND coord_trend = 1, and vice versa.

t_f=mu/kb;

%other params don't change between runs from dataset
t_n = 1e-2; %timescale for neural activity
c_ma = 5; c_prop = 1;  %musc. activity feedback strength, prop feedback strength
a = 1; I = 0; %neural voltage model param, AVB input bias current
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (sech(x-2)).^2; %derivative of sigma(x)
gridsz = 1; %no. of gridpoints per segment
dim = 2; %chain of 2 units
dt = 1e-3;

% ----  I. FIND PERIODIC ORBIT  ----
tic
[ X, period ] = single_oscillator_LC( dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );
toc
fprintf('finding prc')
tic
% ----  II.  CALCULATE iPRC ---- 
Z = single_oscillator_PRC(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime);
elapsedTime = toc
if elapsedTime > 100
    fprintf('PRC took too long to compute, quitting');
    stable_state = NaN;
    yes_it_agrees = 0;
    return
else
% ---- III.  CALCULATE G-FUNCTION  ----
%compute coupling functions - in particular the mechanical coupling fn
[ H_m, ~, ~ ] = oscillator_coupling_fns( X, Z ,dt,period);
theta = 0:.0001:1;
g = -(H_m(-theta)-H_m(theta));

%search g1 for zeros (phase locked states)
pl_states = theta(abs(g)/max(g)<=1e-3);
%find stable phase locked state
pl_states = uniquetol(pl_states,1e-2);
[~, I2] = intersect(theta,pl_states,'stable');
for mm = 1:size(pl_states,2)
    %check if decreasing through phase locked state
    try
    if g(I2(mm)+10)<g(I2(mm))
            stable_state = pl_states(mm);
    end
    catch
       if g(I2(mm))<g(I2(mm)-10)
            stable_state = pl_states(mm);
       end
    end 
end  

%check if stable state is antiphase AND no_trend = 0 or other condition
if abs(stable_state-0.5)<0.05 && ~no_trend
    yes_it_agrees = 1;
elseif abs(stable_state-0.5)<0.05 && no_trend
    yes_it_agrees = 0;
elseif abs(stable_state-0.5)>0.05 && no_trend
    yes_it_agrees = 1;
else
    yes_it_agrees = 0;
end
end
    
end

