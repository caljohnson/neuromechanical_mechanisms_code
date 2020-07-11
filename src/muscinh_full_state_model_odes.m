function [ oderhs ] = muscinh_full_state_model_odes(dim, gridsz, c_ma, c_prop, t_f, ...
    t_n, t_m, a, I, sigma, Gamma, eps_prop, eps_gap, par_no )
%MUSCINH_FULL_STATE_MODEL_ODES  Generates ODEs for the full nervous system 
%                         and body simulation with muscle inhibition
%   Input:   dim - number of oscillators
%            gridsz - no. of body points per oscillator segment
%            c_ma - muscle strength
%            c_prop - local proprioceptive feedback strength
%            t_f - mechanical force timescale
%            t_m - neural activity timescale
%            t_m - muscle activity timescale
%            a - neural parameter
%            I - AVB input bias current
%            sigma - sigma(x) is the thresholding function for the muscle
%                       tension from muscle activity variable
%           Coupling parameters:
%            Gamma - fluid drag coefficient (for mech. coupling), nondiml'zd
%            eps_prop   - strength of nonlocal proprioception
%            eps_gap - strength of electrical NN coupling
%            par_no - oscillator no. to inhibit muscles
%
%  Output: oderhs - ODE RHS functions for the full model


%Proprioceptive Coupling Matrix W (K -> V)
% W = [ c   0 ...
%       e*c c  0 ...
%       0   e*c c 0 ...
%           ...     ...
%                  0  e*c c ]
% nonlocal proprioceptive info flows anterior->posterior
if size(c_prop,2)<=1 %if just one c_prop given, set c_prop to be homogeneous
    c_prop = repmat(c_prop, [1,dim]); %down the oscillator chain
end
if dim>1
    W = zeros(dim,dim);
    W(1,1) = c_prop(1); %local term is full strength
    for ii=2:dim-1
        W(ii, ii) = c_prop(ii); %local term is full strength
        W(ii,ii-1) = -c_prop(ii)*eps_prop; %nonlocal term is eps_prop percent of local term and negative
    end
    W(dim,dim-1) = -c_prop(dim)*eps_prop;%nonlocal term is eps_prop percent of local term & negative
    W(dim,dim) = c_prop(dim);%local term is full strength
else %if dim=1, W is just the local prop. term
    W = c_prop(1);
end

%Electrical coupling matrix E (V -> V)
e = ones(gridsz*dim,1); 
E = eps_gap.*spdiags([e 0*e e], [-1 0 1], dim, dim);


%Mechanical coupling matrix RHS_matrix (K' -> K)
delX = 1/(gridsz*dim); %grid spacing
%2nd difference operator w/ torque-free BCs
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%identity
ey = speye(gridsz*dim,gridsz*dim);
%mechanical PDE (gamma/kb I + mu/kb AA')k_t = -AA'(k + cM)
%Gamma = gamma/kb, tauf = mu/kb
RHS_matrix = (Gamma*ey+t_f*(A*A'))\(-(A*A'));
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
K_dot = @(t,K, AV, AD,c) RHS_matrix*(K + c_ma.*(sigma(S*AV)-sigma(S*AD)));

AV_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_D - volt_V - A);

volt_V_dot = @(t,Kappa,volt_V) (1/t_n).*(volt_V - a.*volt_V.^3 + I +... 
                W*Pa*Kappa + E*volt_V);
volt_D_dot = @(t,Kappa,volt_D) (1/t_n).*(volt_D - a.*volt_D.^3 + I -...
                W*Pa*Kappa + E*volt_D);
            
%paralyze one oscillator by inhibiting muscle activity:
MIer = speye(dim, dim);
if par_no >0
    MIer(par_no, par_no) = 0;
end

%compile all RHS ode functions into one        
% X(1:gridsz*dim) = gridsz*dim K variables
% X(gridsz*dim+1: gridsz*dim+dim) = dim AV variables
% X(gridsz*dim+1+dim: gridsz*dim+2*dim) = dim AD variables
% X(gridsz*dim+1+2*dim: gridsz*dim+3*dim) = dim volt_V variables
% X(gridsz*dim+1+3*dim: gridsz*dim+4*dim) = dim volt_D variables
oderhs = @(t,X) [K_dot(t,X(1:gridsz*dim),X(gridsz*dim+1: gridsz*dim+dim),...
            X(gridsz*dim+1+dim: gridsz*dim+2*dim)); ...
    AV_dot(t,MIer*X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
    MIer*X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1: gridsz*dim+dim));...
     AD_dot(t,MIer*X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
    MIer*X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1+dim: gridsz*dim+2*dim));...
    volt_V_dot(t, X(1:gridsz*dim), X(gridsz*dim+1+2*dim: gridsz*dim+3*dim));...
    volt_D_dot(t, X(1:gridsz*dim), X(gridsz*dim+1+3*dim: gridsz*dim+4*dim));];



end

