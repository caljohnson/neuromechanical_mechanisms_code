% ==================================================================
%                full_model_eps_gamma_loop_zero_muf.m
%
%                      ---------------------------------- 
%   Loop over eps_gap, eps_prop in mu_f=0
%       
% ==================================================================

addpath('./src');
clear

load('fullmodel_eps_gamma_loop.mat')
mu_f = 0;
CN = (3.4*1e-9)*mu_f;

plstate = zeros(size(eps_gaps,2),dim-1)
for kk = 1:size(eps_gaps,2)
 
%make init cond
dt = t_n/10; 

eps_gap = eps_gaps(kk)
eps_prop = eps_props(kk)

%make init cond
phis = 0.7*ones(dim-1,1);
[ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
toc

%make oderhs
oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
    t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );

%run model
TF = 1e2;
%run for TF time, if not enough cycles for phase differences, run again
notlongenough=true;
no_takes = 1; %number of TF runs
while notlongenough
    tic;
    if no_takes==1
    [t,y1] = ode23(oderhs,[0,TF], init_cond);
    else
        %use end of last y as start
        [t2, y2] = ode23(oderhs,[TF*(no_takes-1),TF*no_takes], y1(end,:));
        %stack em together
        t = [t(1:end-1);t2;];
        y1 = [y1(1:end-1,:); y2;];
    end
    toc;
  
    tic;
    %interpolate and compute phase differences
    dt = 1e-3;
    t0 = 0:dt:TF*no_takes;
    y = interp1(t,y1,t0);
    
    try
        %compute phase differences relative to osc 1
        [phis_full,new_period]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
        toc;
         plstate(kk,:) = reshape(mean(phis_full(:,end-10:end),2)',[1 dim-1 1]);
       
        notlongenough = false;
    catch
        no_takes=no_takes+1;
        if no_takes >= 10
            fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Stopping simulation. ']);
            break
        else
            fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
        end
    end
end
end


load('colorblind_colormap.mat');


% fig2 = figure(2); clf;
% colors = colorblind([1 2 6 7 8],:,:);
% h = plot(eps_gaps, plstate,'o-','LineWidth',2,'MarkerSize',10); hold on;
% set(h, {'color'}, num2cell(colors,2));
% ylabel({['pairwise'] ['phase difference']});
% ylim([0.5 1]); xlabel('gap-junctional coupling strength \epsilon_g');
% text(0.51,mod(plstate(end,1),1)-.02,'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
%     text(0.51, mod(plstate(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
%     text(0.51, mod(plstate(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
%     text(0.51,mod(plstate(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
%     text(0.51,mod(plstate(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
% set(gca,'FontSize',30); title('\mu_f = 0');
  