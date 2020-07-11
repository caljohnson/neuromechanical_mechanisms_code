function [Z2] = single_oscillator_adjoint_Eulerstep(Z,X_LC, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime)
%SINGLE_OSCILLATOR_ADJOINT_EULERSTEP Takes an Euler method step for the
% adjoint equations for the single-oscillator model
%
% Inputs:
%       Z - vector of adjoint state values at current time step 
%          Z = (Z_Kappa, Z_A1, Z_A2, Z_V1, Z_V2)
%       X_LC - vector of limit cycle state values at current time step 
%          X_LC = (Kappa, A1, A2, V1, V2)
%       dt - time step size for Euler method step
%       c_ma - muscle strength
%       c_prop - proprioceptive feedback strength
%       t_f - mechanical force timescale
%       t_m - neural activity timescale
%       t_m - muscle activity timescale
%       a - neural parameter
%       sigma_prime - sigma_prime(x) is the derivative of the thresholding 
%               function sigma(x)
%               
%
% Output:
%      Z2 - computed vector of adjoint state values at next time step
%           Z2 = (Z_Kappa, Z_A1, Z_A2, Z_V1, Z_V2)

  %store Jacobian of ODE system
  J = zeros(5,5);
  %Partial derivatives of Kappa eqn:
  J(1,:) = (1/t_f).*[-1,-c_ma.*sigma_prime(X_LC(2)), c_ma*sigma_prime(X_LC(3)), 0, 0];
  %Partial derivatives of A1 eqn:
  J(2,:) = (1/t_m).*[0, -1, 0, 1, -1];
  %Partial derivatives of A2 eqn:
  J(3,:) = (1/t_m).*[0, 0, -1, -1, 1];
  %Partial derivatives of V1 eqn:
  J(4,:) = (1/t_n).*[c_prop, 0, 0, 1-3*a*X_LC(4)^2, 0];
  %Partial derivatives of V2 eqn:
  J(5,:) = (1/t_n).*[-c_prop, 0, 0, 0, 1-3*a*X_LC(5)^2];
  
  %Euler step for adjoint equations
  dZ = -J'*Z;
  Z2=Z+dZ(1:5)*(-dt);

end

