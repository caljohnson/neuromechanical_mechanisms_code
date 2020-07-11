function [ Z ] = single_oscillator_PRC(X_LC, dt, c_ma, c_prop, t_f, t_n, t_m, a, sigma_prime )
%SINGLE_OSCILLATOR_PRC Computes the infinitessimal Phase Response Curves for the
% single-oscillator model
%
% Inputs:
%       X_LC - matrix of limit cycle state values over one period
%          X_LC(t) = (Kappa(t), A1(t), A2(t), V1(t), V2(t))
%       dt - time step size of LC indeces
%       c_ma - muscle strength
%       c_prop - proprioceptive feedback strength
%       t_f - mechanical force timescale
%       t_m - neural activity timescale
%       t_m - muscle activity timescale
%       a - neural parameter
%       sigma_prime - sigma_prime(x) is the derivative of the thresholding 
%               function sigma(x)
%
% Output:
%      Z - computed vector of adjoint state values at each time step in the limit
%             cycle, Z(t)= (Z_Kappa(t), Z_A1(t), Z_A2(t), Z_V1(t), Z_V2(t))
%      period - period of LC (in time units)

ii = size(X_LC,1); %get index-length of limit cycle
y0= ones(5,1); %initial condition 
dy = 10000;
dycrit = 1e-4;
Z = zeros(size(X_LC));

time0 = tic;
timeLimit = 30*1; % 1 minute == 60 seconds
while (dy>dycrit)
% iterate until periodic orbit found to desired precision (dxcrit).
% or if it takes too long
   if toc(time0)>timeLimit
      fprintf('iPRC takes too long'); break
    end
    mm=ii;
    y = y0;
    Z(mm,1:5)=y;            
    while (mm>1)
        mm=mm-1;      
        y=single_oscillator_adjoint_Eulerstep(y,X_LC(mm,1:5), dt, c_ma, c_prop, ...
            t_f, t_n, t_m, a, sigma_prime);
        Z(mm,1:5)=y;
    end;   
    dy=norm(y-y0);
    y0=y;
end;

% normalize Z, so that dXLC/dt*Z=1 for all t. 
dv=diff(X_LC(1:ii,:))/dt;
z0=(Z(1:ii-1,:)+Z(2:ii,:))/2;
sc = zeros(ii-1,1);
for k=1:ii-1
    sc(k)=z0(k,1:5)*transpose(dv(k,1:5));
end;
sc0=median(sc);
Z=Z/sc0;


end

