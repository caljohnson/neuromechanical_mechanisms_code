function [ X ] = phasediffs_to_full_timetrace( dim, phis, X_LC, period, dt, gridsz )
%PHASEDIFFS_TO_FULL_TIMETRACE Creates an time-trace matrix
%   for the full state model from a vector of phase-locked phase differences
%   between each oscillator module
%
% Inputs:
%       dim - number of oscillators in chain
%       phis - vector of phase-locked phase differences
%       X_LC - matrix of single oscillator states over one period
%       period - period of limit cycle in time units
%       dt - time spacing of X_LC vectors
%       gridsz - no. of grid points per oscillator segment 
%      
% Outputs:
%       X - time-trace matrix in full state-variables
%           = [K_1(t); ....; K_dim(t); A_V1(t); ... ; A_Vdim(t); A_D1(t);
%                   ...; A_Ddim(t); V_V1(t); ...; V_Vdim(t); V_D1(t);
%                   ...; V_Ddim(t);];

if size(phis,2)>size(phis,1) %make sure phis is a column vector
    phis = phis';
end
phis = mod(phis,1); %make sure phase goes from 0 to 0.9999
thetas = mod([0; cumsum(phis);],1); %turn phase differences into phases

X = zeros(4*dim+gridsz*dim,10*period/dt); %allocate space
for jj=1:4
    for kk = 1:dim
    %For oscillator kk, set state varitable jj+1 according to the given 
    %limit cycle X_LC at the given initial phase thetas(kk)
    X((jj-1)*dim+kk,:) = repmat([X_LC(round(thetas(kk)*period/dt)+1:end,jj+1);...
                            X_LC(1:round(thetas(kk)*period/dt),jj+1);],[10,1]);
    end
end

%spreading operator - to spread curvature of oscillator to curvatures on 
% the oscillator body segment, keeping a fixed value

% npdf_ids = (-1:2/(gridsz-1):1)'; %for Gaussian Spreading

for kk = 1:dim
     %For oscillator kk, set state varitable 1 according to the given 
    %limit cycle X_LC at the given initial phase thetas(kk)

    X(4*dim+(kk-1)*gridsz+1:4*dim+kk*gridsz,:) = repmat(ones(gridsz,1)*[X_LC(round(thetas(kk)*period/dt)+1:end,1);...
                            X_LC(1:round(thetas(kk)*period/dt),1);]',[1,10]);
 %with Gaussian spreading   
% X(4*dim+(kk-1)*gridsz+1:4*dim+kk*gridsz,:) = repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids))*[X_LC(round(thetas(kk)*period/dt)+1:end,1);...
%                             X_LC(1:round(thetas(kk)*period/dt),1);]',[1,10]);
end

%rearrange so curvatures are first
X = [X(4*dim+1:end,:); X(1:4*dim,:);];

end

