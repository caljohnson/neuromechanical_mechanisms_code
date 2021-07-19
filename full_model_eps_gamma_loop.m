% ==================================================================
%                full_model_eps_gamma_loop.m
%
%                      ---------------------------------- 
%   Loop over eps_gap  
%   then find eps_prop to get correct wavelength in water
%   then loop over external viscosities and examine coordination
%       
% ==================================================================

addpath('./src');
clear

% % -- NM MODEL PARAMETERS --
% mu = 1.3e-7;
% kb = 2.6e-7;
% t_f=mu/kb;
% t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
% c_ma = 10; c_prop = 1;  %musc. activity feedback strength, prop feedback strength
% a = 1; I = 0; %neural voltage model param, AVB input bias current
% nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature
% sigma = @(x) (1/2).*tanh(x-2); %thresholding function for muscle activity into tension
% sigma_prime = @(x) (1/2).*(sech(x-2)).^2; %derivative of sigma(x)
% gridsz = 1; %no. of gridpoints per segment
% dim = 6; 
% %external fluid viscosities
% mu_f = [1, 10, logspace(2,3,4), 10^4, 10^4.4472];
% 
% %vary eps_gap
% % eps_gaps = [0.01:.01:0.05 0.1:.1:.5];
% eps_gaps = [0 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1];
% 
% periods = zeros(size(mu_f,2),size(eps_gaps,2));
% wavelns = zeros(size(mu_f,2),size(eps_gaps,2));
% mean_amp = zeros(size(mu_f,2),size(eps_gaps,2));
% plstate = zeros(size(mu_f,2),dim-1,size(eps_gaps,2));
% eps_props = zeros(size(eps_gaps));
% 
% for kk = 1:size(eps_gaps,2)
%  
% %make init cond
% dt = t_n/10; 
% 
% eps_gap = eps_gaps(kk)
% %loop to determine t_m so that the period is within .1 of 0.5 sec
% period = 0;
% true_period = 0.5;
% t_m_low = 0.05;
% t_m_hi = 0.25;
% no_runs = 1; %number of times run this loop
% %bisection on t_m to find in [50 ms ,250 ms]
% while abs(period - true_period)>10^-1 && no_runs < 20
%     tic
%     t_m = (t_m_hi+t_m_low)/2
%     [ X_LC, period ] = single_oscillator_LC( dt, c_ma(1), c_prop, t_f, t_n, t_m, a, I, sigma );
%     period
%     no_runs = no_runs+1;
%     if period < true_period
%         t_m_low = t_m;
%     else
%         t_m_hi = t_m;
%     end
% end
% %make init cond
% phis = 0.7*ones(dim-1,1);
% [ init_cond ] = phases_to_init_cond( dim, phis, X_LC, period, dt, gridsz );
% toc
% 
% %first loop to determine eps_prop and eps_gap so that the low-viscosity
% %wavelength is within 0.1 of 1.5 wavelengths/bodylength
% jj=1;
% CN = (3.4*1e-9)*mu_f(jj);
% 
% wvln = 0;
% true_wvln = 1.54;
% % eps_prop = 0.05;
% % eps_step = .01;
% % no_runs = 1;
% % while abs(wvln - true_wvln)>5*10^-2 && no_runs < 20
% %     if no_runs>1
% %         if wvln>true_wvln
% %             eps_prop = eps_prop+eps_step
% %             eps_step = eps_step*.9
% %         else
% %             eps_prop = abs(eps_prop-eps_step)
% %             eps_step = eps_step*.9
% %         end
% %     end
% eps_prop = eps_gap*2.7;
% eps_prop_low = max(1e-5, eps_gap*2.7e-1);
% eps_prop_hi = min(eps_gap*2.7e2,1);
% no_runs = 1;
% while abs(wvln - true_wvln)>5*10^-2 && no_runs < 15
%     if no_runs>1
%         if wvln>true_wvln
%             eps_prop_low = eps_prop;
%             
%         else
%             eps_prop_hi = eps_prop;
%         end
%         eps_prop= (eps_prop_hi+eps_prop_low)/2
%     end
% 
%     tic
%     %make oderhs
%     oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
%                 t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
%     %run model
%     TF = 3e2;
%     %run for TF time, if not enough cycles for phase differences, run again
%     notlongenough=true;
%     no_takes = 1; %number of TF runs
%     while notlongenough
%         tic;
%         if no_takes==1
%             [t,y1] = ode23(oderhs,[0,TF], init_cond);
%         else
%             %use end of last y as start
%             [t2, y2] = ode23(oderhs,[TF*(no_takes-1),TF*no_takes], y1(end,:));
%             %stack em together
%             t = [t(1:end-1);t2;];
%             y1 = [y1(1:end-1,:); y2;];
%         end
%         toc;
%   
%         tic;
%         %interpolate and compute phase differences
%         dt = 1e-3;
%         t0 = 0:dt:TF*no_takes;
%         y = interp1(t,y1,t0);
%     
%         try
%             %compute phase differences relative to osc 1
%             [phis_full,new_period]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
%             toc;
%             plstate(jj,:,kk) = reshape(mean(phis_full(:,end-10:end),2)',[1 dim-1 1]);
%             periods(jj,kk) = new_period*dt;
%             freq = 1./periods(jj);
%             wavelns(jj,kk) = 1/((6/5)*sum(1-mod(plstate(jj,:,kk),1)));
%             wvln = wavelns(jj,kk)
%             notlongenough = false;
%         catch
%             no_takes=no_takes+1;
%             fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
%         end
%     end
%     no_runs = no_runs+1;
% end
% 
% eps_props(kk) = eps_prop
% if no_runs<20
% 
% %loop over rest of gammas to show gait adaptation
% for jj = 2:size(mu_f,2)
%     CN = (3.4*1e-9)*mu_f(jj);
% %make oderhs
% oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
%     t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
% 
% %run model
% TF = 1e2;
% %run for TF time, if not enough cycles for phase differences, run again
% notlongenough=true;
% no_takes = 1; %number of TF runs
% while notlongenough
%     tic;
%     if no_takes==1
%     [t,y1] = ode23(oderhs,[0,TF], init_cond);
%     else
%         %use end of last y as start
%         [t2, y2] = ode23(oderhs,[TF*(no_takes-1),TF*no_takes], y1(end,:));
%         %stack em together
%         t = [t(1:end-1);t2;];
%         y1 = [y1(1:end-1,:); y2;];
%     end
%     toc;
%   
%     tic;
%     %interpolate and compute phase differences
%     dt = 1e-3;
%     t0 = 0:dt:TF*no_takes;
%     y = interp1(t,y1,t0);
%     
%     try
%         %compute phase differences relative to osc 1
%         [phis_full,new_period]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
%         toc;
%          plstate(jj,:,kk) = reshape(mean(phis_full(:,end-10:end),2)',[1 dim-1 1]);
%          periods(jj,kk) = new_period*dt;
%          wavelns(jj,kk) = 1/((6/5)*sum(1-mod(plstate(jj,:,kk),1)));
% 
%         notlongenough = false;
%     catch
%         no_takes=no_takes+1;
%         if no_takes >= 10
%             fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Stopping simulation. ']);
%             break
%         else
%             fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
%         end
%     end
% end
% end
% 
% else
%     disp(strcat('wavelength in water not achievable for eps_gap',num2str(eps_gap)));
% end
% end
% 
% fprintf('finished all trials');
% legendCell = cellstr(strcat('\',num2str(eps_gaps', 'epsilon_g=%.5f')));
% legendCell{end+1} = 'FY 2010';
% 
% 
% %-- Loop over eps_gap and eps_prop, find wavelength in water
% jj=1;
% CN = (3.4*1e-9)*mu_f(jj);
% eps_props = 2.7*eps_gaps;
% wavelns2 = zeros(size(eps_props,2),size(eps_gaps,2));
% for kk=1:size(eps_gaps,2)
%     eps_gap = eps_gaps(kk)
% for jj = 1:size(eps_props,2)
%     tic
%     eps_prop = eps_props(jj)
%     %make oderhs
%     oderhs = full_state_model_odes(dim, gridsz, c_ma, c_prop, mu,kb, ...
%                 t_n, t_m, a, I, sigma, CN, eps_prop, eps_gap );
%     %run model
%     TF = 1e2;
%     %run for TF time, if not enough cycles for phase differences, run again
%     notlongenough=true;
%     no_takes = 1; %number of TF runs
%     while notlongenough
%         tic;
%         if no_takes==1
%             [t,y1] = ode23(oderhs,[0,TF], init_cond);
%         else
%             %use end of last y as start
%             [t2, y2] = ode23(oderhs,[TF*(no_takes-1),TF*no_takes], y1(end,:));
%             %stack em together
%             t = [t(1:end-1);t2;];
%             y1 = [y1(1:end-1,:); y2;];
%         end
%         toc;
%   
%         tic;
%         %interpolate and compute phase differences
%         dt = 1e-3;
%         t0 = 0:dt:TF*no_takes;
%         y = interp1(t,y1,t0);
%     
%         try
%             %compute phase differences relative to osc 1
%             [phis_full,new_period]  = full_timetrace_to_phasediffs( dim, y, gridsz,20 );
%             toc;
%             plstate = mean(phis_full(:,end-10:end),2)';
%             wavelns2(jj,kk) = 1/((6/5)*sum(1-mod(plstate,1)));
%             notlongenough = false;
%         catch
%             no_takes=no_takes+1;
%             fprintf(['TF = ' num2str(TF*no_takes) ' not long enough. Running for longer. ']);
%         end
%     end
% end
% end
% save('fullmodel_eps_gamma_loop_new.mat')

load('fullmodel_eps_gamma_loop_new.mat')
% load('colorblind_colormap.mat');
% colors = colorblind([1 2 4:8 10 12],:,:);
colors = cool(9);
figure(1); clf;
h = semilogx(mu_f,wavelns(:,1:4),'LineWidth',4.0,'Marker','+','Markersize',10); hold on;
set(h, {'color'}, num2cell(colors(1:4,:),2));
h2 = semilogx(mu_f,wavelns(:,5:end),'LineWidth',4.0,'Marker','o','Markersize',10); hold on;
set(h2, {'color'}, num2cell(colors(5:end,:),2));
semilogx(mu_f, [1.54, 1.375, 1.25, 1.2, 1, .9, .8, .75],':k','LineWidth',3,'Marker','x','Markersize',20);
set(gca,'FontSize',30);
legend(legendCell,'location','northeastoutside','FontSize',20);
xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
ylabel({['wavelength'] ['(\lambda / L)']});

% 
% figure(2); clf;
% loglog(eps_gaps,eps_props,'LineWidth',4.0,'Marker','o','Markersize',10); hold on;
% p = polyfit(eps_gaps(1:5),eps_props(1:5),1); 
% f = polyval(p,eps_gaps(1:9)); 
% loglog(eps_gaps(1:9),f,':','LineWidth',3.0); 
% xlabel('\epsilon_g'); ylabel('\epsilon_p');
% legend('fit to wavelength in water', strcat('linear, m = ',num2str(p(1))));
% set(gca,'FontSize',30);

% for jj=1:10
%     fig1 = figure(jj+2); clf;
%     colors = colorblind([1 2 6 7 8],:,:);
%     h = semilogx(mu_f,mod(plstate(:,:,jj),1),'o-','LineWidth',2,'MarkerSize',10); hold on;
%     set(h, {'color'}, num2cell(colors,2));
%     ylabel({['pairwise'] ['phase difference']});
%     ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
%     set(gca,'FontSize',30); 
%     title(['\epsilon_g = ', num2str(eps_gaps(jj)), ', \epsilon_p = ', num2str(eps_props(jj))])
%     text(3.5*10^4,mod(plstate(end,1,jj),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
%     text(3.5*10^4, mod(plstate(end,2,jj),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
%     text(3.5*10^4, mod(plstate(end,3,jj),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
%     text(3.5*10^4,mod(plstate(end,4,jj),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
%     text(3.5*10^4,mod(plstate(end,5,jj),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
%     saveas(fig1,['figs/eps_gap_loop_figs/eps_gap_phasediffs_' num2str(jj) '.png']);
% end




% save('fullmodel_eps_gamma_loop_new.mat')
% eps_props = repmat(eps_props,size(eps_props,2),1)
% n_prop = size(eps_props,2);
% for jj=1:size(eps_props,1)
%    eps_props(jj,n_prop+1) = eps_gaps(jj)*2.7; 
% end

figure(2); clf;
h = semilogx(eps_props,wavelns2(:,1:4),'LineWidth',4.0,'Marker','+','Markersize',10); hold on;
set(h, {'color'}, num2cell(colors(1:4,:),2));
h2 = semilogx(eps_props,wavelns2(:,5:end),'LineWidth',4.0,'Marker','o','Markersize',10); hold on;
set(h2, {'color'}, num2cell(colors(5:end,:),2));
set(gca,'FontSize',30);
legend(legendCell,'location','northeastoutside','FontSize',20);
ylabel({['wavelength'] ['(\lambda / L)']});
xlabel('proprioceptive coupling strength \epsilon_p');
set(gca, 'FontSize',30); ylim([0,2]); xlim([2e-4, 10^0]);
