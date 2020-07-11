function [  ] = timetrace_to_curv_kymo( Xt, dim, gridsz, fig_no )
%TIMETRACE_TO_CURV_KYMO Creates a curvature kymograph
%   from the time-trace of the full-state model
%
%   Inputs:
%       Xt - matrix of model states over time
%       dim - number of oscillators in chain
%       gridsz - no. of grid points per oscillator segment 
%       fig_no - figure number to plot on

%transpose Xt if need be
if size(Xt,1)>size(Xt,2)
    Xt = Xt';
end

Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
figure(fig_no); clf;
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
colormap(blueblackred); caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
ylabel('time'); xlabel('body position'); 
set(gca,'FontSize',30)
end

