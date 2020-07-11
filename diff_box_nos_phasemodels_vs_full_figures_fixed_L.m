% ==================================================================
%
%                     diff_box_nos_phasemodels_vs_full_figures_fixed_L.m
%                      ------- 
%  Runs the full model and phase model for the N-box chain of NM
%  oscillators to investigate coupling
%       with fixed body length L


% and generates the figures of N-oscillator wavelengths and phase-diffs
% ==================================================================
addpath('./src');
addpath('./Nbox_phasediffs/fixed_L');
load('colorblind_colormap.mat')

%2 boxes
load('2box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/2box_phasediffs.png');

%3 boxes
load('3box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/3box_phasediffs.png');
% text(3.5*10^4,0.775,'\phi_3', 'FontSize',30,'Color',colors(3,:,:));
% text(3.5*10^4,0.625,'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));

%turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(2); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths'] ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);



%4 boxes
load('4box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/4box_phasediffs.png');
% text(3.5*10^4,0.625,'\phi_4', 'FontSize',30,'Color',colors(4,:,:));
% text(3.5*10^4,0.825,'\phi_5', 'FontSize',30,'Color',colors(5,:,:));

% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(4); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);


%5 boxes
load('5box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/5box_phasediffs.png');

% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(6); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%6 boxes
load('sixbox_osc_vs_full_gamma_loop_data.mat');
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
text(3.5*10^4,mod(phase_model_eq2(end,5),1),'\phi_4', 'FontSize',30,'Color',colors(5,:,:));
saveas(fig,'figs/Nbox_phasediffs/6box_phasediffs.png');

% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(6); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%7 boxes
load('7box_phase_vs_full_gamma_loop_data.mat');
fig=figure(7); clf;
colors = colorblind([1 2 4 6 7 8 ],:,:);
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
text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
saveas(fig,'figs/Nbox_phasediffs/7box_phasediffs.png');

% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(8); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%8 boxes
load('8box_phase_vs_full_gamma_loop_data.mat');
fig=figure(8); clf;
colors = colorblind([1 2 4 6 7 8 9],:,:);
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
text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
saveas(fig,'figs/Nbox_phasediffs/8box_phasediffs.png');

% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(10); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%9 boxes
load('9box_phase_vs_full_gamma_loop_data.mat');
fig=figure(9); clf;
colors = colorblind([1 2 4 6 7 8 9 10],:,:);
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
text(3.5*10^4,mod(phase_model_eq2(end,6),1),'\phi_6', 'FontSize',30,'Color',colors(6,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,7),1),'\phi_7', 'FontSize',30,'Color',colors(7,:,:));
text(3.5*10^4,mod(phase_model_eq2(end,8),1),'\phi_8', 'FontSize',30,'Color',colors(8,:,:));
saveas(fig,'figs/Nbox_phasediffs/9box_phasediffs.png');


% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(12); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%10 boxes
load('10box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/10box_phasediffs.png');

% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%11 boxes
load('11box_phase_vs_full_gamma_loop_data.mat');
fig=figure(11); clf;
% colors = colorblind([1 2 3 4 6 7 8 9 10 12],:,:);
colors = linspecer(10);
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
text(3.5*10^4,mod(phase_model_eq2(end,10),1),'\phi_{10}', 'FontSize',30,'Color',colors(10,:,:));
saveas(fig,'figs/Nbox_phasediffs/11box_phasediffs.png');
% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%12 boxes
load('12box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/12box_phasediffs.png');
% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%15 boxes
load('15box_phase_vs_full_gamma_loop_data.mat');
fig=figure(15); clf;
colors = linspecer(14);
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
text(3.5*10^4,mod(phase_model_eq2(end,14),1),'\phi_{14}', 'FontSize',30,'Color',colors(14,:,:));
saveas(fig,'figs/Nbox_phasediffs/15box_phasediffs.png');
% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%20 boxes
load('20box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/20box_phasediffs.png');
% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);

%25 boxes
load('25box_phase_vs_full_gamma_loop_data.mat');
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
saveas(fig,'figs/Nbox_phasediffs/25box_phasediffs.png');
% 
% %turn into wavelength
% wvln2 = 1./((6/(dim-1))*sum(1-mod(phase_model_eq2,1),2));
% figure(14); clf; str1=['\epsilon_{g} = ' num2str(eps_gap) ' and \epsilon_{p} = ' num2str(eps_prop)];
% semilogx(little_gamma,wvln2,'o-','LineWidth',2,'MarkerSize',10); hold on;
% wvlns_full2 = 1./((6/(dim-1))*sum(1-mod(plstate2,1),2));
% semilogx(little_gamma,wvlns_full2,'+','LineWidth',2,'MarkerSize',10);
% xlabel('External Fluid Viscosity \mu_f (mPa \cdot s)'); 
% ylabel({['wavelengths']  ['(\lambda/L)']});
% set(gca,'FontSize',30); title(str1);