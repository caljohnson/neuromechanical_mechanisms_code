% ==================================================================
%                make_mukb_trend_colormap_figs.m
%
%  Uses data generated in paperfig_parametrize_tmtf_loop_regimes.m to make figs
%                      ---------------------------------- 
%  Figure of good and bad kymographs
%   Figure of tau_b vs. mu_b grid of coordination outcomes
%   Figures of example good/bad wavelength trends
%       
% ==================================================================


% load('trials_data.mat')
addpath('./src');
% load('paper_trials_data.mat');
% load('paper_trials_data_new.mat');
load('paper_trials_data_4.mat');

%initialize empty vectors
N = size(trial,2);

kbs = zeros(N,6);
mus = zeros(N,6);
tau_ms = zeros(N,6);
trends = zeros(N,6);
tau_fs = zeros(N,6);
periods = zeros(N,6);
wvlns = zeros(N,6,8);
pdiffs = zeros(N,6,8,5);
% twobox_agreement = ones(N,6);
% twobox_plstate = zeros(N,6,8);
twobox_phaseagreement = ones(N,6);
twobox_phasestate = zeros(N,6);
% sixbox_agreement = ones(N,6);
% sixbox_pdiffs = zeros(N,6,8,5);
% sixbox_phase_wvlns = zeros(N,6,8);

for jj=1:N
for kk = 1:6
    kbs(jj,kk) = trial(jj).run(kk).k_b;
    mus(jj,kk) = trial(jj).run(kk).mu;
    tau_ms(jj,kk) = trial(jj).run(kk).tau_m;
    trends(jj,kk) = trial(jj).run(kk).no_trend;
    tau_fs(jj,kk) = trial(jj).run(kk).tau_f;
    periods(jj,kk) = trial(jj).run(kk).single_osc_period;
    wvlns(jj,kk,:) = trial(jj).run(kk).wavelengths;
    pdiffs(jj,kk,:,:) = trial(jj).run(kk).phase_diffs;
    twobox_phaseagreement(jj,kk) = trial(jj).run(kk).phase_agrees;
    twobox_phasestate(jj,kk) = trial(jj).run(kk).phase_state; 
    
    %decide whether coordination is coherent (0) or not (1)
    %check wavelength trend
    coord_trends(jj,kk) = 0;
    not_traveling(jj,kk) = 0;
    eps_wvln = 1;
    for mm=1:size(wvlns(jj,kk,:),3)-1
        if wvlns(jj,kk,mm)<wvlns(jj,kk,mm+1)-eps_wvln
            coord_trends(jj,kk) = 1;
        break
        end
    end
    %also check phase differences for spread - indicating no traveling wave
    for mm=1:size(wvlns(jj,kk,:),3)-1
    if max(pdiffs(jj,kk,mm,:)) - min(pdiffs(jj,kk,mm,:)) >= 0.5
%         coord_trends(jj,kk) = 1;
        not_traveling(jj,kk) = 1;
        break
    end
    end
    
        
end
end

I2 = find(not_traveling);
I1 = setdiff(find(coord_trends), I2);
I3 = setdiff(find(~not_traveling), I1);

J = find(abs(periods-0.5)>0.1);


M = find(abs(twobox_phasestate - 0.5) > 0.1);


figure(1); clf; colormap(gray(1))
scatter(log10(mus(I2)),log10(tau_fs(I2)),250,'r^','filled'); hold on
scatter(log10(mus(I1)),log10(tau_fs(I1)),250,'m+','LineWidth',4); hold on
% scatter(log10(mus(I1_v2)),log10(tau_fs(I1_v2)),250,'r^','filled'); hold on
scatter(log10(mus(I3)),log10(tau_fs(I3)),250,'bo', 'filled');
scatter(log10(mus(J)),log10(tau_fs(J)),400,'ks','filled');

% scatter(log10(mus(nottravel_inds)),log10(tau_fs(nottravel_inds)),400,'g+');
% scatter(log10(mus(M)),log10(tau_fs(M)),400,'m*');
% contour(log10(mus),log10(tau_fs),(tau_ms),'--','LineWidth', '4.0');
%  [c,h] = contour(log10(mus),log10(tau_fs),(tau_ms),[0.0594 0.0875 0.125 0.125 0.1625 0.1813],'--');
%  h.LineWidth = 3;
xlabel('log \mu_b'); ylabel('log \tau_b');
legend('no traveling wave', 'incorrect wavelength trend', 'correct wavelength trend', ...
    'incorrect frequency', 'Location', 'northoutside');%, 'phase model predicts synchrony');

 set(gca,'FontSize',30)
%  clabel(c,h, 'FontSize',20, 'LineWidth',4)
% scatter(log10(mus(K)),log10(tau_fs(K)),500,'kx');
% scatter(log10(mus(L)),log10(tau_fs(L)),300,'ro');
%, '\tau_m contours');
% legend('incorrect wavelength trend', 'correct wavelength trend', ...
%     'incorrect frequency','2-osc. phase model disagrees', '6-osc. phase model disagrees');


%-----------------------------------------
%box a few cases to display in next figure

%just the 2 good, 1 bad
circle_mus = [log10(mus(end-6,4));  log10(mus(end-9,2)); log10(mus(1,end));];
circle_tfs = [log10(tau_fs(end-6,4));  log10(tau_fs(end-9,2)); log10(tau_fs(1,end));];

% for full grid
% circle_mus = [log10(mus(end-6,2)); log10(mus(end-6,4)); log10(mus(end-6,6));...
%     log10(mus(end-9,2)); log10(mus(end-9,4)); log10(mus(end-9,6));...
%     log10(mus(end-12,2)); log10(mus(end-12,4)); log10(mus(end-12,6));];
% circle_tfs = [log10(tau_fs(end-6,2)); log10(tau_fs(end-6,4)); log10(tau_fs(end-6,6));...
%     log10(tau_fs(end-9,2)); log10(tau_fs(end-9,4)); log10(tau_fs(end-9,6));...
%     log10(tau_fs(end-12,2)); log10(tau_fs(end-12,4)); log10(tau_fs(end-12,6));];
scatter(circle_mus, circle_tfs, 400, 'ks');


%taum vs tau
%-------------- tau_m labels
[c,h] = contour(log10(mus),log10(tau_fs),(tau_ms),[0.0563 0.1 0.15 0.2 0.24],':k');
 h.LineWidth = 3;
 xlim([min(min(log10(mus))) max(max(log10(mus)))]);
% text(-6.8,0.5,{'Incorrect', 'Frequency'}, 'FontSize',20);
 text(-6.8,-0.07143,'\tau_m = .05 s', 'FontSize',20);
 text(-6.8,-0.25,'\tau_m = .10 s', 'FontSize',20);
 text(-6.8,-.7,'\tau_m = .15 s', 'FontSize',20);
%  text(-6.75,-1.2,'\tau_b \approx \tau_m = .15 s', 'FontSize',30,'Rotation',90);
 text(-6.8,-1.15,'\tau_m = .20 s', 'FontSize',20);
 text(-6.8,-1.6,'\tau_m = .24 s', 'FontSize',20);
 set(gca,'FontSize',30)
 set(gcf,'Position',[1    59   640   646]);
 
 %------- tau vs tau
  figure(10); clf; 

  loglog(tau_ms(I2), tau_fs(I2),'r^', 'MarkerSize',15,'MarkerFaceColor','r'); hold on;
loglog(tau_ms(I3), tau_fs(I3),'bo', 'MarkerSize',20,'MarkerFaceColor','b');
loglog(tau_ms(I1), tau_fs(I1),'m+', 'MarkerSize',10,'MarkerFaceColor','m');
loglog(tau_ms(J), tau_fs(J),'ks','MarkerSize',20,'MarkerFaceColor','k');
legend('no traveling wave',...
    'correct wavelength trend', ...
    'incorrect wavelength trend', 'incorrect frequency');%,...
%     '\tau_f = \tau_m', '\tau_f = 1');
loglog(logspace(-2,0,10),logspace(-2,0,10),'k--', 'LineWidth',4.0);
loglog(logspace(-2,0,10),ones(10,1),'k--', 'LineWidth',4.0);

xlabel('log \tau_m'); ylabel('log \tau_b');
set(gca,'FontSize',30) 
 
  
 %--------- show kymograph of traveling wave and not traveling wave

little_gamma = trial(1).run(1).ext_gamma;
fig3 = figure(3); clf;
t_n = 1e-2; %timescales for neural (10ms)
c_ma = 5;%musc. activity feedback strength
c_prop = 1;  %prop feedback strength (arbitrary)
a = 1; I = 0; %neural voltage model param, AVB input bias current
sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
gridsz = 1; dim = 6; %body discretization parameters
dt0 = 0.1*t_n;
 
 %traveling wave (high \gamma)
 [X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, tau_fs(end-6,6), t_n, tau_ms(end-6,6), a, I, sigma );
Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(6).phase_diffs(end,:), X_LC, period, dt0, gridsz );
subplot(1,2,1); 
Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([round(size(Kappa,2)/2) size(Kappa,2)]);
colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
ylabel('time \rightarrow'); 
xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20,...
    'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
title({'Traveling wave', strcat('(log \tau_b = ',num2str(log10(tau_fs(end-6,6)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-6,6)),'%.2f'), ', \mu_f = 28 Pa s)')}); %set(gca,'FontSize',20)
    
    
%not traveling wave kymo
%gamma high
% end-6,2 TOP LEFT in fig2
[X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, tau_fs(end-6,2), t_n, tau_ms(end-6,2), a, I, sigma );
Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(2).phase_diffs(end,:), X_LC, period, dt0, gridsz );
subplot(1,2,2); 
Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([round(size(Kappa,2)/2) size(Kappa,2)]);
colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
ylabel('time \rightarrow'); 
set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
xlabel('body coordinate');
title({'No traveling wave', strcat('(log \tau_b = ',num2str(log10(tau_fs(end-6,2)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-6,2)),'%.2f'), ', \mu_f = 28 Pa s)')}); 

%---------------- plot 2 good wavelength trends, 1 bad one -------
fig2 = figure(2); clf; 
% subplot(1,3,1);
semilogx(little_gamma, reshape(wvlns(end-6,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor','b', 'Color','b'); 
ylim([0,1.6]); ylabel({['wavelength'] ['(\lambda / L)']});
%ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
%   title(strcat('\tau_b = ',num2str((tau_fs(end-6,4)),'%.2f'),...
%         ' s , \mu_b = ', num2str((mus(end-6,4)),'%1.1e'),'N(mm^2)s'));
set(gca,'FontSize',30)
% title('(a)');
  xlabel('\mu_f (mPa \cdot s)'); 
  
fig4 = figure(4); clf;
semilogx(little_gamma, reshape(wvlns(end-9,2,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor', 'b', 'Color','b');
    ylim([0,1.6]); ylabel({['wavelength'] ['(\lambda / L)']});
%     title('(b)');
%       title(strcat('\tau_b = ',num2str((tau_fs(end-9,2)),'%.2f'),...
%         ' s , \mu_b = ', num2str((mus(end-9,2)),'%1.1e'), ' N(mm^2)s'));
set(gca,'FontSize',30)
  xlabel('\mu_f (mPa \cdot s)'); 
    
fig5 = figure(5);
semilogx(little_gamma, reshape(wvlns(1,end,:),[1,8]),':','LineWidth',4,'Marker','+','Markersize',10,'MarkerFaceColor','m', 'Color','m'); 
ylim([0,3]);  ylabel({['wavelength'] ['(\lambda / L)']});
xlabel('\mu_f (mPa \cdot s)');  %title('(c)');
%   title(strcat('\tau_b = ',num2str((tau_fs(1,end)),'%.2f'),...
%         ' s, \mu_b = ', num2str((mus(1,end)),'%1.1e'), ' N(mm^2)s'));
set(gca,'FontSize',30)

return
    
% 
% 
%  %-------------------------------------------%-------------------------------------------
%  %grid of sample wavelength trends
%  little_gamma = trial(1).run(1).ext_gamma;
% fig2 = figure(2); clf;
% subplot(3,3,1); semilogx(little_gamma, reshape(wvlns(end-6,2,:),[1,8]),':','LineWidth',4,'Marker','^','Markersize',10,'MarkerFaceColor','r', 'Color','r');
%     ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
%     title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,2)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,2)),'%.2f'))); set(gca,'FontSize',20)
%     
% subplot(3,3,2); semilogx(little_gamma, reshape(wvlns(end-6,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor','b', 'Color','b'); 
% ylim([0,1.6]);
%   title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,4)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,4)),'%.2f')));set(gca,'FontSize',20)
% 
% subplot(3,3,3); semilogx(little_gamma, reshape(wvlns(end-6,6,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor','b', 'Color','b'); 
% ylim([0,1.6]);
%   title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,6)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,6)),'%.2f')));set(gca,'FontSize',20)
% 
% subplot(3,3,4); semilogx(little_gamma, reshape(wvlns(end-9,2,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor', 'b', 'Color','b');
%     ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
%       title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,2)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-9,2)),'%.2f')));set(gca,'FontSize',20)
% 
% subplot(3,3,5); semilogx(little_gamma, reshape(wvlns(end-9,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor', 'b', 'Color','b'); 
% ylim([0,1.6]);
%  title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,4)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-9,4)),'%.2f')));set(gca,'FontSize',20)
%     
% subplot(3,3,6); semilogx(little_gamma, reshape(wvlns(end-9,6,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor','b','Color','b'); 
% ylim([0,1.6]);
%  title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,4)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-9,4)),'%.2f')));set(gca,'FontSize',20)
%     
% subplot(3,3,7); semilogx(little_gamma, reshape(wvlns(end-12,2,:),[1,8]),':','LineWidth',4,'Marker','^','Markersize',10,'MarkerFaceColor','r', 'Color','r'); 
%     ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
%     xlabel('\mu_f (mPa \cdot s)'); 
%      title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,2)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-12,2)),'%.2f')));set(gca,'FontSize',20)
%     
% subplot(3,3,8); semilogx(little_gamma, reshape(wvlns(end-12,4,:),[1,8]),':','LineWidth',4,'Marker','^','Markersize',10,'MarkerFaceColor','r', 'Color','r'); 
% ylim([0,1.6]);
% xlabel('\mu_f (mPa \cdot s)'); 
%      title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,4)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-12,4)),'%.2f')));set(gca,'FontSize',20)
%     
% subplot(3,3,9); semilogx(little_gamma, reshape(wvlns(end-12,6,:),[1,8]),':','LineWidth',4,'Marker','^','Markersize',10,'MarkerFaceColor','r', 'Color','r'); 
% ylim([0,1.6]);
%      title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,6)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-12,6)),'%.2f')));set(gca,'FontSize',20)
% % semilogx(little_gamma, [1.54, 1.375, 1.25, 1.2, 1, .9, .8, .75],'o');
% xlabel('\mu_f (mPa \cdot s)'); 
% set(gcf,'Position',[1          59        1280         646]);
% 
% %-------------------------------------------%-------------------------------------------
% %grid of kymographs
% %grid of sample wavelength trends
% little_gamma = trial(1).run(1).ext_gamma;
% fig3 = figure(3); clf;
% t_n = 1e-2; %timescales for neural (10ms)
% c_ma = 5;%musc. activity feedback strength
% c_prop = 1;  %prop feedback strength (arbitrary)
% a = 1; I = 0; %neural voltage model param, AVB input bias current
% sigma = @(x) tanh(x-2); %thresholding function for muscle activity into tension
% gridsz = 1; dim = 6; %body discretization parameters
% dt0 = 0.1*t_n;
% 
% % end-6,2 TOP LEFT in fig2
% [X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, tau_fs(end-6,2), t_n, tau_ms(end-6,2), a, I, sigma );
% % gamma low
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(2).phase_diffs(1,:), X_LC, period, dt0, gridsz );
% subplot(3,3,1); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); caxis([min(Kappa(1,:)) max(Kappa(1,:))]); %colorbar();
% % ylabel('time'); xlabel('body position'); 
% title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,2)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,2)),'%.2f'))); set(gca,'FontSize',20)
%     ylabel('\mu_f = 1 mPa s'); set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment', 'right');
%     set(gca,'YTickLabel',[], 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% %gamma med
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(2).phase_diffs(4,:), X_LC, period, dt0, gridsz );
% subplot(3,3,4); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
% ylabel('\mu_f = 215 mPa s'); set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment', 'right');
% set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% % ylabel('time'); xlabel('body position'); 
% %gamma high
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(2).phase_diffs(end,:), X_LC, period, dt0, gridsz );
% subplot(3,3,7); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
% % ylabel('time'); 
% ylabel('\mu_f = 28 Pa s'); xlabel('body position'); 
% set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment', 'right');
% set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% 
% 
% % end-6,4 top middle in fig2
% [X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, tau_fs(end-6,4), t_n, tau_ms(end-6,4), a, I, sigma );
% % gamma low
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(4).phase_diffs(1,:), X_LC, period, dt0, gridsz );
% subplot(3,3,2); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
% % ylabel('time'); xlabel('body position'); 
% title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,4)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,4)),'%.2f')));
%    set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]); 
% %gamma med
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(4).phase_diffs(4,:), X_LC, period, dt0, gridsz );
% subplot(3,3,5); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
% set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% % ylabel('time'); xlabel('body position'); 
% %gamma high
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(4).phase_diffs(end,:), X_LC, period, dt0, gridsz );
% subplot(3,3,8); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); %colorbar();
% % ylabel('time'); 
% xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% 
% 
% % end-6,6 top right in fig2
% [X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, tau_fs(end-6,6), t_n, tau_ms(end-6,6), a, I, sigma );
% % gamma low
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(6).phase_diffs(1,:), X_LC, period, dt0, gridsz );
% subplot(3,3,3); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
% % ylabel('time'); xlabel('body position'); 
% title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,6)),'%.2f'),...
%         ', log \mu_b = ', num2str(log10(mus(end-6,6)),'%.2f')));
%     set(gca,'YTickLabel',[],'FontSize',20,  'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% %gamma med
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(6).phase_diffs(4,:), X_LC, period, dt0, gridsz );
% subplot(3,3,6); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
% % ylabel('time');
% % xlabel('body position'); 
% set(gca,'YTickLabel',[],'FontSize',20, 'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% %gamma high
% Xt = phasediffs_to_full_timetrace( dim, 1-trial(end-6).run(6).phase_diffs(end,:), X_LC, period, dt0, gridsz );
% subplot(3,3,9); 
% Kappa = Xt(1:gridsz*dim,:);
% Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
% surf(Kappa');
% view(2); shading flat;
% xlim([1 gridsz*dim+1]); ylim([0, size(Kappa,2)]);
% colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
% % ylabel('time'); 
% xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20,...
%     'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
% % colorbar(

