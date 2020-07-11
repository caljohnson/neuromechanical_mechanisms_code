function [ LC, period ] = single_oscillator_LC( dt0, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma )
%SINGLE_OSCILLATOR_LC Computes the limit cycle for the 
% single-oscillator model
%
% Inputs:
%       dt0 - time step size for limit cycle vector
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
%
% Output:
%      LC - computed vector of state values at each time step in the limit
%             cycle, LC(t)= (Kappa(t), A1(t), A2(t), V1(t), V2(t))
%      period - period of LC (in time units)

nv = 5; %number of state variables
nmax=1e8;  % maximal number of time steps in one period
v = zeros(nmax,nv);
Kmark=0.01;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt=dt0;
x0(1:nv)=[1;0;0;1;-1]; % initial state (K,A1,A2,V1,V2)

dx=10000.0;     
dxcrit=1e-5; % precision for periodic orbit

time0 = tic;
timeLimit = 10*1; % in seconds (1 minute == 60 seconds)

while (dx>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
% or if it takes too long
   if toc(time0)>timeLimit
      fprintf('finding LC takes too long'); 
      LC = NaN;
%       LC=v(1:ii,:);
      period = NaN;
      return 
    end
    x=x0;
    ii=1;             
    flagg=0;
    while (flagg==0 && ii<=nmax) % time-step until system returns to
                                 % v=vmark with dv/dt>0
        y=single_oscillator_Eulerstep(x,dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma); 
        if toc(time0)>timeLimit
         fprintf('finding LC takes too long'); break
        end
        if (x(1)<Kmark && y(1)>=Kmark && ii>1)  
            dt=(Kmark-x(1))/(y(1)-x(1))*dt;
            y=single_oscillator_Eulerstep(x,dt, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma);      
            flagg=1; 
            dt=dt0;
        end;
        x=y;
        v(ii,1:nv)=x(1:nv); %output
        ii=ii+1;
    end;    
    dx=norm(x-x0);
    x0=x;
end;
ii=ii-1; %set ii to actual period index-length
period = ii*dt; %period in time-units

% rearrange output so that max v is at t=0.
[~,imax]=max(v(1:ii,1));
v2(1:ii-imax+1,1:nv)=v(imax:ii,1:nv);
v2(ii-imax+2:ii,1:nv)=v(1:imax-1,1:nv);
LC=v2;

end

