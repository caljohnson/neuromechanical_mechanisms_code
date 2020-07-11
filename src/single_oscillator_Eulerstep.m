function [ Y ] = single_oscillator_Eulerstep(X, dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma)
%SINGLE_OSCILLATOR_EULERSTEP Takes an Euler method step for the
% single-oscillator model
%
% Inputs:
%       X - vector of state values at current time step 
%          X = (Kappa, A1, A2, V1, V2)
%       dt - time step size for Euler method step
%       c_ma - muscle strength
%       c_prop - proprioceptive feedback strength
%       t_f - mechanical force timescale
%       t_m - neural activity timescale
%       t_m - muscle activity timescale
%       a - neural parameter
%       I - AVB input bias current
%       sigma - sigma(x) is the thresholding function for the muscle
%               tension from muscle activity variable
%
% Output:
%      Y - computed vector of state values at next time step
%           Y = (Kappa, A1, A2, V1, V2)

  % -- K --
  Y(1) = X(1) + (-X(1) - c_ma.*(sigma(X(2)) - sigma(X(3))))*dt/t_f;  
   
  % -- A1 --
  Y(2) = X(2) + (-X(2) + X(4) - X(5))*dt/t_m;
  
  % -- A2 --
  Y(3) = X(3) + (-X(3) + X(5) - X(4))*dt/t_m;
  
  % -- V1 --
  Y(4) = X(4) + (X(4)-a*X(4)^3 + I + c_prop*X(1))*dt/t_n;
  
  % -- V2 --
  Y(5) = X(5) + (X(5)-a*X(5)^3 + I - c_prop*X(1))*dt/t_n;

     
end

