function [ pl_states, state_eigs, stabilities ] = phase_locked_state_solve( dim,phase_odes)
%PHASE_LOCKED_STATE_SOLVE Solves for the phase-locked states Phi* of the
%       oscillator chain using the phase model ODES
%           Phi' = F(Phi) -> F(Phi*) = 0
%       and their stabilities via the eigenvalues of the jacobian
%
% Inputs:
%       dim - number of oscillators in chain
%       phase_odes - RHS ode functions for the phase model
%           Phi' = F(Phi)
%      
% Outputs:
%       pl_states - matrix of unique phase-locked states
%       stabilities - matrix of corresponding stabilities

phi0 = (0:.1:1)';

foptions = optimoptions('fsolve','Display', 'none', 'MaxFunctionEvaluations', 1e3);
tol = 2e-1;
%set empty arrays to store states and stabilities in
pl_states = [];

if dim == 4
    %loop over initial phi0 in 3-dim space
for j1 = 1:size(phi0,1)
    for j2 = 1:size(phi0,1)
        for j3 = 1:size(phi0,1)
           %set initial condition for fsolve
           x0 = [phi0(j1);  phi0(j2); phi0(j3);];         
           x = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
           %check if phase-locked state is unique up to tol
           pl_states = uniquetol([pl_states; mod(x',1);], tol, 'ByRows', 1);
        end
    end
    %go through unique states, check & store stability
    stabilities = zeros(size(pl_states,1),3);
    state_eigs = zeros(size(pl_states,1),3);
        for jj = 1:size(pl_states,1)
            x0 = pl_states(jj,:)';
            [~,~,~,~,J] = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
            state_eigs(jj,:) = eig(J);
            stabilities(jj,:) = (state_eigs(jj,1)<0)*(state_eigs(jj,2)<0)*(state_eigs(jj,3)<0);
        end
end
elseif dim==3
    for j1 = 1:size(phi0,1)
        for j2 = 1:size(phi0,1)
           %set initial condition for fsolve
           x0 = [phi0(j1);  phi0(j2);];         
           x = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
           %check if phase-locked state is unique up to tol
           pl_states = uniquetol([pl_states; mod(x',1);], tol, 'ByRows', 1);
        end
    end
    %go through unique states, check & store stability
    stabilities = zeros(size(pl_states,1),2);
    state_eigs = zeros(size(pl_states,1),2);
        for jj = 1:size(pl_states,1)
            x0 = pl_states(jj,:)';
            [~,~,~,~,J] = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
            state_eigs(jj,:) = eig(J);
            stabilities(jj,:) = (state_eigs(jj,1)<0)*(state_eigs(jj,2)<0);
        end

elseif dim==2
        for j1 = 1:size(phi0,1)
           %set initial condition for fsolve
           x0 = phi0(j1);
           x = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
           %check if phase-locked state is unique up to tol
           pl_states = uniquetol([pl_states; mod(x,1);], tol, 'ByRows', 1);
        end
        %go through unique states, check & store stability
        stabilities = zeros(size(pl_states,1),1);
        state_eigs = zeros(size(pl_states,1),1);
        for jj = 1:size(pl_states,1)
            x0 = pl_states(jj);
            [~,~,~,~,J] = fsolve(@(phi) phase_odes(1,phi), x0,foptions);
            state_eigs(jj) = eig(J);
            stabilities(jj) = (state_eigs(jj)<0);
        end
end


    

end

