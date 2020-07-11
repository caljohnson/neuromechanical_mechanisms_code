function [ phis,period_cycle ] = full_timetrace_to_phasediffs( dim, X, gridsz,N)
%FULL_TIMETRACE_TO_PHASEDIFFS Creates an matrix of phase-differences
%between each oscillator module in the full-state model over time
%
% Inputs:
%       dim - number of oscillators in chain
%       X - time-trace of full-state model
%       gridsz - no. of grid points per oscillator segment 
%       N - number of periods to moving-window over
%      
% Outputs:
%       phis - time-trace matrix in phase-diffs
%           = [phi_1(t); ....; phi_dim-1(t);];

if size(X,2)<size(X,1) %transpose y so that state-variables are rows
    X = X';
end

%compute average curvature of each oscillator - Pa is averaging matrix
rowinds = repelem(1:dim,gridsz);
colinds = 1:gridsz*dim;
npdf_ids = -1:2/(gridsz-1):1;
% Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));
Pa = sparse(rowinds, colinds, (1/gridsz).*repmat(normpdf(npdf_ids)/max(normpdf(npdf_ids)),1,dim));
Kappa = Pa*X(1:gridsz*dim,:);

%flip Kappa so that each time series is a row 
if size(Kappa,2) < size(Kappa,1)
    Kappa = Kappa';
end

%get autocorrelation of first oscillator
control_corr =  xcorr(Kappa(1,:));
%find peaks of autocorr
[vals,locs] = findpeaks(control_corr);
[~, sorted_locs] = sort(vals,'descend');
%peak-to-peak distance defines period of oscillation
period_cycle = abs(locs(sorted_locs(1))-locs(sorted_locs(2)));

%use moving window of size N periods
window = N*period_cycle;
slide_count = round(size(Kappa,2)/period_cycle);
phis = zeros(dim-1,slide_count);

for tt = 1:slide_count
    
    %compute phase differences in this window
    inds = (tt-1)*period_cycle+1:(tt-1)*period_cycle+window;
    %stop when window reaches end of timeseries
    if inds(end) > size(Kappa,2)
        phis(:,tt:end) = [];
        break
    end
    
    for ii=2:dim
        %get autocorrelation of first oscillator
         control_corr =  xcorr(Kappa(ii-1,inds));
         %find peaks of autocorr
        [vals,locs] = findpeaks(control_corr);
        [~, sorted_locs] = sort(vals,'descend');
        %peak-to-peak distance defines period of oscillation
        period_cycle = abs(locs(sorted_locs(1))-locs(sorted_locs(2)));
        %get max-peak index of autocorr
        [~,I_control] = max(control_corr);
        %compute cross-correlation between oscillator ii and ii-1 
        cross_cor = xcorr(Kappa(ii-1,inds), Kappa(ii,inds));
        %get index of max-peak of cross-corr
        [~,I_new] = max(cross_cor);
        %phase difference is peak-peak distance / cycle period
        phis(ii-1,tt) = (I_new-I_control)/period_cycle;
   
        %ensure phase-diff is between 0 and 1
        phis(ii-1,tt) = mod(phis(ii-1,tt),1);
       
    end
end


end

