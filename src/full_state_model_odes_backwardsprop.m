function [ oderhs ] = full_state_model_odes_backwardsprop(dim, gridsz, c_ma, c_prop, mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap )
%FULL_STATE_MODEL_ODES  Generates ODEs for the full nervous system 
%                         and body simulation with backwards proprioception
%   Input:   dim - number of oscillators
%            gridsz - no. of body points per oscillator segment
%            c_ma - muscle strength
%            c_prop - local proprioceptive feedback strength
%            mu - internal viscosity
%            kb - internal elasticity
%            t_m - neural activity timescale
%            t_m - muscle activity timescale
%            a - neural parameter
%            I - AVB input bias current
%            sigma - sigma(x) is the thresholding function for the muscle
%                       tension from muscle activity variable
%           Coupling parameters:
%            CN - fluid drag coefficient
%            eps_prop   - strength of nonlocal proprioception
%            eps_gap - strength of electrical NN coupling
%
%  Output: oderhs - ODE RHS functions for the full model


%Proprioceptive Coupling Matrix W (K -> V)
% W = [ c   e ...
%       0 c  e ...
%       0   0 c e ...
%           ...     ...
%                  0  0 c ]
% nonlocal proprioceptive info flows backwards posterior->anterior
if size(c_prop,2)<=1 %if just one c_prop given, set c_prop to be homogeneous
    c_prop = repmat(c_prop, [1,dim]); %down the oscillator chain
end
if dim>1
    W = zeros(dim,dim);
    W(1,1) = c_prop(1); %local term is full strength
    for ii=1:dim-1
        W(ii, ii) = c_prop(ii); %local term is full strength, positive
        W(ii,ii+1) = eps_prop; %nonlocal, posterior term is weak, positive
    end
    W(dim,dim) = c_prop(dim);%local term is full strength, positive
else %if dim=1, W is just the local prop. term
    W = c_prop(1);
end

%Electrical coupling matrix E (V -> V)
e = ones(gridsz*dim,1); 
E = spdiags([e -2*e e], [-1 0 1], dim, dim);
E(1,1) = -1; E(end,end) = -1; %no left/right neighbors at head/tail
E = eps_gap.*E;

%Mechanical coupling matrix RHS_matrix (K' -> K)
L0 = 1; %initial length of whole worm (in mm), 1mm=1e-3 m
delX = L0/(gridsz*dim); %grid spacing
%2nd difference operator
D2 = spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%4th difference operator
D4 = D2*D2';
D4(1,1) = 7; D4(end,end) = 7; %moment-free, force-free BCs
D4 = (1/delX^4).*D4;
%identity
ey = speye(gridsz*dim,gridsz*dim);
%mechanical PDE (CN I + mu d4)k_t = - kb D4(k + cM)
RHS_matrix = (CN*ey+mu*D4)\(-kb*D4);
RHS_matrix(abs(RHS_matrix)<1e-16) = 0;

%spread operator - for spreading muscle tension to the body segments
rowinds = 1:gridsz*dim;
colinds = repelem(1:dim, gridsz);
npdf_ids = -1:2/(gridsz-1):1;
S = sparse(rowinds, colinds, repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));
%proprioceptive kernel
rowinds = repelem(1:dim,gridsz);
colinds = 1:gridsz*dim;
% Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));
Pa = sparse(rowinds, colinds, (1/gridsz).*repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));

%ODEs
K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_ma.*(sigma(S*AV)-sigma(S*AD)));

AV_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_D - volt_V - A);

volt_V_dot = @(t,Kappa,volt_V) (1/t_n).*(volt_V - a.*volt_V.^3 + I +... 
                W*Pa*Kappa + E*volt_V);
volt_D_dot = @(t,Kappa,volt_D) (1/t_n).*(volt_D - a.*volt_D.^3 + I -...
                W*Pa*Kappa + E*volt_D);

%compile all RHS ode functions into one        
% X(1:gridsz*dim) = gridsz*dim K variables
% X(gridsz*dim+1: gridsz*dim+dim) = dim AV variables
% X(gridsz*dim+1+dim: gridsz*dim+2*dim) = dim AD variables
% X(gridsz*dim+1+2*dim: gridsz*dim+3*dim) = dim volt_V variables
% X(gridsz*dim+1+3*dim: gridsz*dim+4*dim) = dim volt_D variables
oderhs = @(t,X) [K_dot(t,X(1:gridsz*dim),X(gridsz*dim+1: gridsz*dim+dim),...
            X(gridsz*dim+1+dim: gridsz*dim+2*dim)); ...
    AV_dot(t,X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
    X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1: gridsz*dim+dim));...
     AD_dot(t,X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
    X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1+dim: gridsz*dim+2*dim));...
    volt_V_dot(t, X(1:gridsz*dim), X(gridsz*dim+1+2*dim: gridsz*dim+3*dim));...
    volt_D_dot(t, X(1:gridsz*dim), X(gridsz*dim+1+3*dim: gridsz*dim+4*dim));];



end

