function [ H_m, H_p, H_g ] = oscillator_coupling_fns( X, Z ,dt,period)
%OSCILLATOR_COUPLING_FNS Coupled oscillator model coupling functions:
%       Creates anonymous functions of x that take any input x
%       and mod it by 1, that defines the coupling in terms of phase
%       from 0 to 1 for each coupling form
%
% Inputs:
%       X - matrix of limit cycle state values
%          X = (Kappa, A1, A2, V1, V2)
%       Z - matrix of phase response values
%          Z = (Z_Kappa, Z_A1, Z_A2, Z_V1, Z_V2)
%       dt - time step size for Euler method step
%       period - period length in time units
%      
% Outputs:
%       H_m - coupling function for mechanical coupling (K' -> K)
%       H_p - coupling function for proprioceptive coupling (K -> V)
%       H_g - coupling function for electrical coupling (V -> V)


%K' from limit cycle
Kdv = (-[X(3:end,1); X(1:2,1)]+...
    8*[X(2:end,1); X(1,1)]-8*[X(end,1); X(1:end-1,1)]+...
    [X(end-1:end,1); X(1:end-2,1)])/(12*dt);

%---- calculate H-functions --------
h_m_phi = -ifft(-fft(flip(Z(:,1))).*fft(Kdv))*dt/period; %mechanical coupling (all-to-all)
h_p_phi = -ifft(-fft(flip(Z(:,4))).*fft(X(:,1)))*dt/period; %proprioceptive (one-way) coupling
h_g_phi = -ifft(-fft(flip(Z(:,4))).*fft(X(:,4)))*dt/period;%neural (gap-junction) coupling

%turn H-functions into MATLAB functions
tees = linspace(0,1,size(h_m_phi,1))';

H_m = fit(tees, h_m_phi,'smoothingspline');%,'SmoothingParam',0.99);
H_m = @(t) H_m(mod(t,1));

H_p = fit(tees, h_p_phi,'smoothingspline');%,'SmoothingParam',0.99);
H_p = @(t) H_p(mod(t,1));

H_g = fit(tees, h_g_phi,'smoothingspline');%,'SmoothingParam',0.99);
H_g = @(t) H_g(mod(t,1));


end

