% ==================================================================
%
%                     diff_box_nos_phasemodels_vs_full_figures_fixed_ell.m
%                      ------- 
%  From data of the full model and phase model for the N-box chain of NM
%       with fixed module length ell


%  generates the figures of N-oscillator wavelengths and phase-diffs
% ==================================================================
addpath('./src');
load('colorblind_colormap.mat')

%2 boxes
load('Nbox_phasediffs/fixed_ell/2box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(2); clf;
colors = colorblind([1],:,:);
% plot phase model states
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
% plot full model states
h2 = semilogx(little_gamma,plstate2,'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/2box_phasediffs_v2.png');

mm=1;

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1;

%3 boxes
load('Nbox_phasediffs/fixed_ell/3box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(3); clf;
colors = colorblind([1 2],:,:);
% plot phase model states
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
% plot full model states
h2 = semilogx(little_gamma,plstate2,'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/3box_phasediffs_v2.png');
% text(3.5*10^4,0.775,'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,0.625,'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%4 boxes
load('Nbox_phasediffs/fixed_ell/4box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(4); clf;
colors = colorblind([1 2 6],:,:);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/4box_phasediffs_v2.png');
% text(3.5*10^4,0.625,'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1


%5 boxes
load('Nbox_phasediffs/fixed_ell/5box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(5); clf;
colors = colorblind([1 2 6 7],:,:);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/5box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%6 boxes
load('Nbox_phasediffs/fixed_ell/6box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(6); clf;
colors = colorblind([1 2 6 7 8],:,:);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/6box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%10 boxes
load('Nbox_phasediffs/fixed_ell/10box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(10); clf;
% colors = colorblind([1 2 4 6 7 8 9 10 12],:,:);
colors = linspecer(9);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,8),1),'\phi_8', 'FontSize',30,'Color',colors(8,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,9),1),'\phi_9', 'FontSize',30,'Color',colors(9,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/10box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%12 boxes
load('Nbox_phasediffs/fixed_ell/12box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(12); clf;
% colors = colorblind([1 2 3 4 6 7 8 9 10 11 12],:,:);
colors = linspecer(11);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,8),1),'\phi_8', 'FontSize',30,'Color',colors(8,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,9),1),'\phi_9', 'FontSize',30,'Color',colors(9,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,10),1),'\phi_{10}', 'FontSize',30,'Color',colors(10,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,11),1),'\phi_{11}', 'FontSize',30,'Color',colors(11,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/12box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%20 boxes
load('Nbox_phasediffs/fixed_ell/20box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(20); clf;
colors = linspecer(19);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,8),1),'\phi_8', 'FontSize',30,'Color',colors(8,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,9),1),'\phi_9', 'FontSize',30,'Color',colors(9,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,10),1),'\phi_{10}', 'FontSize',30,'Color',colors(10,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,11),1),'\phi_11', 'FontSize',30,'Color',colors(11,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,12),1),'\phi_12', 'FontSize',30,'Color',colors(12,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,13),1),'\phi_13', 'FontSize',30,'Color',colors(13,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,14),1),'\phi_14', 'FontSize',30,'Color',colors(14,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,15),1),'\phi_{15}', 'FontSize',30,'Color',colors(15,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,16),1),'\phi_16', 'FontSize',30,'Color',colors(16,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,17),1),'\phi_17', 'FontSize',30,'Color',colors(17,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,18),1),'\phi_18', 'FontSize',30,'Color',colors(18,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,19),1),'\phi_{19}', 'FontSize',30,'Color',colors(19,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/20box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));
mm= mm+1

%25 boxes
load('Nbox_phasediffs/fixed_ell/25box_phase_vs_full_gamma_loop_data_v2.mat');
fig=figure(25); clf;
colors = linspecer(24);
h = semilogx(little_gamma,mod(phase_model_eq2,1),'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,mod(plstate2,1),'+','LineWidth',2,'MarkerSize',10);
set(h2, {'color'}, num2cell(colors,2));
% xlabel('External Fluid Viscosity \gamma (mPa \cdot s)'); 
ylabel({['pairwise'] ['phase difference']});
ylim([0.5 1]); xlabel('External Fluid Viscosity \mu_f mPa s'); 
set(gca,'FontSize',30);
text(3.5*10^4,mod(phase_model_eq2(end,1),1),'\phi_1', 'FontSize',30,'Color',colors(1,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,2),1),'\phi_2', 'FontSize',30,'Color',colors(2,:,:));
% text(3.5*10^4, mod(phase_model_eq2(end,3),1),'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,4),1),'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_5', 'FontSize',30,'Color',colors(5,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,8),1),'\phi_8', 'FontSize',30,'Color',colors(8,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,9),1),'\phi_9', 'FontSize',30,'Color',colors(9,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,10),1),'\phi_{10}', 'FontSize',30,'Color',colors(10,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,11),1),'\phi_11', 'FontSize',30,'Color',colors(11,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,12),1),'\phi_12', 'FontSize',30,'Color',colors(12,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,13),1),'\phi_13', 'FontSize',30,'Color',colors(13,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,14),1),'\phi_14', 'FontSize',30,'Color',colors(14,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,15),1),'\phi_{15}', 'FontSize',30,'Color',colors(15,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,16),1),'\phi_16', 'FontSize',30,'Color',colors(16,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,17),1),'\phi_17', 'FontSize',30,'Color',colors(17,:,:));
% text(3.5*10^4,mod(phase_model_eq2(end,18),1),'\phi_18', 'FontSize',30,'Color',colors(18,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,24),1),'\phi_{24}', 'FontSize',30,'Color',colors(24,:,:));
saveas(fig,'figs/Nbox_phasediffs/fixed_ell/25box_phasediffs_v2.png');

%turn into wavelength
wavelengths_phasemodel(mm,:) = 1/6./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel(mm,:) = 1/6./((sum(1-plstate2,2))./(dim-1));
wavelengths_phasemodel_v2(mm,:) = 1/dim./(sum(1-mod(phase_model_eq2,1),2)./(dim-1));
wavelengths_fullmodel_v2(mm,:) = 1/dim./((sum(1-plstate2,2))./(dim-1));


%wavelength plot of all cases
figure(100); clf; 
colors = colorblind([1 2 4 5 6 7 8 10 12],:,:);
% str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
h = semilogx(little_gamma,wavelengths_phasemodel,'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,wavelengths_fullmodel,'+','LineWidth',2,'MarkerSize',10); 
set(h2, {'color'}, num2cell(colors,2));
xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
ylabel({['wavelengths'] ['(\lambda/L)']});
set(gca,'FontSize',30); %title(str1);
ylim([0,2]);
legend('N=2', 'N=3', 'N=4','N=5','N=6','N=10','N=12','N=20','N=25','location','northeastoutside')

figure(101); clf; 
colors = colorblind([1 2 4 5 6 7 8 10 12],:,:);
% str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
h = semilogx(little_gamma,wavelengths_phasemodel_v2,'o-','LineWidth',2,'MarkerSize',10); hold on;
set(h, {'color'}, num2cell(colors,2));
h2 = semilogx(little_gamma,wavelengths_fullmodel_v2,'+','LineWidth',2,'MarkerSize',10); 
set(h2, {'color'}, num2cell(colors,2));
xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
ylabel({['wavelengths'] ['(\lambda/L_N)']});
set(gca,'FontSize',30); %title(str1);
% ylim([0,6]);
ylim([0,3]);
legend('N=2', 'N=3', 'N=4','N=5','N=6','N=10','N=12','N=20','N=25','location','northeastoutside')



