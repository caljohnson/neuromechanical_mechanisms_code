function [ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz )
%PHASES_TO_INIT_COND Creates an initial condition vector
%   for the full state model from a vector of phase-diffs for the oscillator
%   modules
%
% Inputs:
%       dim - number of oscillators in chain
%       phis - vector of initial phase-differences
%       X_LC - matrix of single oscillator states over one period
%       period - period of limit cycle in time units
%       dt - time spacing of X_LC vectors
%       gridsz - no. of grid points per oscillator segment 
%      
% Outputs:
%       init_cond - initial condition vector in full state-variables
%           = [K_1(0); ....; K_dim(0); A_V1(0); ... ; A_Vdim(0); A_D1(0);
%                   ...; A_Ddim(0); V_V1(0); ...; V_Vdim(0); V_D1(0);
%                   ...; V_Ddim(0);];

if size(phis,2)>size(phis,1) %make sure phis is a column vector
    phis = phis';
end
phis = mod(phis,1); %make sure phase goes from 0 to 0.9999
thetas = mod([0; cumsum(phis);],1); %turn phase differences into phases
init_cond = zeros(4*dim+gridsz*dim,1); %allocate space
for jj=1:4
    %For each oscillator, set state variable jj according to the given 
    %limit cycle X_LC at the given initial phase in the vector thetas
    init_cond((jj-1)*dim+1:jj*dim) = X_LC(mod(round(thetas*period/dt),size(X_LC,1))+1,jj+1)';
end

%spreading operator - to spread curvature of oscillator to a Gaussian of
%curvatures on the oscillator body segment
rowinds = 1:gridsz*dim;
colinds = repelem(1:dim, gridsz);
npdf_ids = -1:2/(gridsz-1):1;
S = sparse(rowinds, colinds, repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));
init_cond(4*dim+1:end) = (S*X_LC(mod(round(thetas*period/dt),size(X_LC,1))+1,1))';

%rearrange so curvatures are first
init_cond = [init_cond(4*dim+1:end); init_cond(1:4*dim);];

end

