
addpath('./src');
load('paper_trials_data_4.mat');
load('colormap_mukb_trends.mat');

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

%Fang-Yen 2010 wavelength data:
FY_wavelengths = interp1(1:8,[1.54, 1.375, 1.25, 1.2, 1, .9, .8, .75]',1:0.25:8);

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
    x = reshape(wvlns(jj,kk,:), [8 1]);
    monotonicity(jj,kk) =  1/7*sum(sign(x(2:end)-x(1:end-1)));
    coord_trends(jj,kk) = 1*(monotonicity(jj,kk)>=-.15);
    match_FY(jj,kk) = immse(interp1(1:8,x,1:.25:8), FY_wavelengths);
    eps_wvln = 0.05;
    
%     for mm=1:size(wvlns(jj,kk,:),3)-3
%         if wvlns(jj,kk,mm)<wvlns(jj,kk,mm+1)-eps_wvln
%             coord_trends(jj,kk) = 1;
%         break
%         end
%     end
    %also check phase differences for spread - indicating no traveling wave
%     for mm=1:size(wvlns(jj,kk,:),3)-3
%         if max(pdiffs(jj,kk,mm,:)) - min(pdiffs(jj,kk,mm,:)) >= 0.5
%             coord_trends(jj,kk) = 1;
%             break
%         end
%     end
end
end

I1 = find(coord_trends);
I2 = find(~coord_trends);
J = find(abs(periods-0.5)>0.1);
% K = find(~twobox_agreement);
% L = find(~sixbox_agreement);
% M1 = find(abs(twobox_phasestate - 0.5) > 0.1);
% M = intersect(M1,~J);
% save('all_loop_data_paper.mat');
% return
fig1= figure(1); clf;
% colormap(colormap(cmap);); caxis([min(log10(match_FY(:))) max(log10(match_FY(:)))]);
colormap(cmap);
scatter(log10(mus(:)),log10(tau_fs(:)),250,log10(match_FY(:)),'filled'); hold on
% scatter(log10(mus(:)),log10(tau_fs(:)),250,monotonicity(:),'filled'); hold on
% scatter(log10(mus(I1)),log10(tau_fs(I1)),250,'ro','filled'); hold on
% scatter(log10(mus(I2)),log10(tau_fs(I2)),250,'b^', 'filled');
scatter(log10(mus(J)),log10(tau_fs(J)),400,'ks','filled');
% cb = colorbar('northoutside'); title(cb,'Monotonicity of Wavelength Trend'); 
cb = colorbar('northoutside'); title(cb,{'Mean-Squared Error of Wavelength Trend' '(vs. Fang-Yen et al. 2010)'}); 
circle_mus = [log10(mus(end-6,1)); log10(mus(end-6,4)); log10(mus(end-6,6));...
    log10(mus(end-9,1)); log10(mus(end-9,4)); log10(mus(end-9,6));...
    log10(mus(end-12,1)); log10(mus(end-12,4)); log10(mus(end-12,6));];
circle_tfs = [log10(tau_fs(end-6,1)); log10(tau_fs(end-6,4)); log10(tau_fs(end-6,6));...
    log10(tau_fs(end-9,1)); log10(tau_fs(end-9,4)); log10(tau_fs(end-9,6));...
    log10(tau_fs(end-12,1)); log10(tau_fs(end-12,4)); log10(tau_fs(end-12,6));];
scatter(circle_mus, circle_tfs, 400, 'ks');

% figure(2); clf;
% scatter(log10(tau_ms(:)),log10(tau_fs(:)),250,monotonicity(:),'filled'); 
% scatter(log10(mus(M)),log10(tau_fs(M)),400,'mx','LineWidth',3);
% contour(log10(mus),log10(tau_fs),(tau_ms),'--','LineWidth',4);
% colormap(gray(1))
[c,h] = contour(log10(mus),log10(tau_fs),(tau_ms),[0.0563 0.1 0.15 0.2 0.24],':k');
 h.LineWidth = 3;
% plot(linspace(-9.7,-6.5,10),-0.8239*ones(10),'k--', 'LineWidth',4.0);
xlabel('log \mu_b (mPa \cdot s)'); ylabel('log \tau_b (s)');
% lgd= legend('monotonicity of wavelength trend', ...
%     'correct frequency in water unattainable');%, 'phase model predicts synchrony');
% lgd.Location = 'northoutside';

xlim([min(min(log10(mus))) max(max(log10(mus)))]);
text(-6.8,0.5,{'Incorrect', 'Frequency'}, 'FontSize',20);
 text(-6.8,-0.07143,'\tau_m = .05 s', 'FontSize',20);
 text(-6.8,-0.25,'\tau_m = .10 s', 'FontSize',20);
 text(-6.8,-.7,'\tau_m = .15 s', 'FontSize',20);
%  text(-6.75,-1.2,'\tau_b \approx \tau_m = .15 s', 'FontSize',30,'Rotation',90);
 text(-6.8,-1.15,'\tau_m = .20 s', 'FontSize',20);
 text(-6.8,-1.6,'\tau_m = .24 s', 'FontSize',20);
 set(gca,'FontSize',30)
 set(gcf,'Position',[1    59   640   646]);
 saveas(fig1,'figs/wvln_trends_tfs_vs_mus.png');
%  clabel(c,h, 'FontSize',20, 'LineWidth',4)
% scatter(log10(mus(K)),log10(tau_fs(K)),500,'kx');
% scatter(log10(mus(L)),log10(tau_fs(L)),300,'ro');
%, '\tau_m contours');
% legend('incorrect wavelength trend', 'correct wavelength trend', ...
%     'incorrect frequency','2-osc. phase model disagrees', '6-osc. phase model disagrees');

% figure(2); clf;
% scatter(log10(tau_ms(:)), log10(tau_fs(:)),100,monotonicity(:),'filled'); hold on
% % loglog(tau_ms(I1), tau_fs(I1),'ro', 'MarkerSize',15,'MarkerFaceColor','r'); 
% % loglog(tau_ms(I2), tau_fs(I2),'b*', 'MarkerSize',10,'MarkerFaceColor','b');
% % loglog(tau_ms(M), tau_fs(M),'mx', 'MarkerSize',20,'MarkerFaceColor','m');
% plot(linspace(-2,0,10),linspace(-2,0,10),'k--', 'LineWidth',4.0);
% % plot(linspace(-2,0,10),zeros(10,1),'k--', 'LineWidth',4.0);
% % legend('monotonicity of wavelength trend', ...
% %     'phase model predicts synchrony',...
% %     '\tau_f = \tau_m', '\tau_f = 1');
% xlabel('log \tau_m'); ylabel('log \tau_f');
% set(gca,'FontSize',30)

little_gamma = trial(1).run(1).ext_gamma;
fig3 = figure(3); clf;
subplot(3,3,1); semilogx(little_gamma, reshape(wvlns(end-6,1,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(17,:), 'Color',cmap(17,:));
    ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
    title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,1)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-6,1)),'%.2f'))); set(gca,'FontSize',20)
    
subplot(3,3,2); semilogx(little_gamma, reshape(wvlns(end-6,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(12,:), 'Color',cmap(12,:)); 
ylim([0,1.6]);
  title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,4)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-6,4)),'%.2f')));set(gca,'FontSize',20)

subplot(3,3,3); semilogx(little_gamma, reshape(wvlns(end-6,6,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(5,:), 'Color',cmap(5,:)); 
ylim([0,1.6]);
  title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-6,6)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-6,6)),'%.2f')));set(gca,'FontSize',20)

subplot(3,3,4); semilogx(little_gamma, reshape(wvlns(end-9,2,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor', cmap(15,:),'Color',cmap(15,:));
    ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
      title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,2)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-9,2)),'%.2f')));set(gca,'FontSize',20)

subplot(3,3,5); semilogx(little_gamma, reshape(wvlns(end-9,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(5,:), 'Color',cmap(5,:)); 
ylim([0,1.6]);
 title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,4)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-9,4)),'%.2f')));set(gca,'FontSize',20)
    
subplot(3,3,6); semilogx(little_gamma, reshape(wvlns(end-9,6,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(1,:),'Color',cmap(1,:)); 
ylim([0,1.6]);
 title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-9,4)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-9,4)),'%.2f')));set(gca,'FontSize',20)
    
subplot(3,3,7); semilogx(little_gamma, reshape(wvlns(end-12,1,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(23,:), 'Color',cmap(23,:)); 
    ylim([0,1.6]); ylabel({['wavelength'] ['per bodylength'] ['(\lambda / L)']});
    xlabel('\mu_f (mPa \cdot s)'); 
     title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,1)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-12,1)),'%.2f')));set(gca,'FontSize',20)
    
subplot(3,3,8); semilogx(little_gamma, reshape(wvlns(end-12,4,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(12,:), 'Color',cmap(12,:)); 
ylim([0,1.6]);
xlabel('\mu_f (mPa \cdot s)'); 
     title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,4)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-12,4)),'%.2f')));set(gca,'FontSize',20)
    
subplot(3,3,9); semilogx(little_gamma, reshape(wvlns(end-12,6,:),[1,8]),':','LineWidth',4,'Marker','o','Markersize',10,'MarkerFaceColor',cmap(12,:), 'Color',cmap(12,:)); 
ylim([0,1.6]);
     title(strcat('log \tau_b = ',num2str(log10(tau_fs(end-12,6)),'%.2f'),...
        ', log \mu_b = ', num2str(log10(mus(end-12,6)),'%.2f')));set(gca,'FontSize',20)
% semilogx(little_gamma, [1.54, 1.375, 1.25, 1.2, 1, .9, .8, .75],'o');
xlabel('\mu_f (mPa \cdot s)'); 
set(gcf,'Position',[1          59        1280         646]);
saveas(fig3,'figs/grid_of_wvln_trends.png');
% ylabel('wavelength per bodylength (\lambda / L)');
% legend('model', 'Fang-Yen et al. 2010');