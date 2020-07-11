% ==================================================================
%
%                     oscillator_properties.m
%
%  This program plots the oscillator cycle properties 
%   period, curvature amplitude as functions of cma and cprop
%
% ==================================================================

addpath('./src');
clear

% %loop over cma, cprop
% c_ma = [0:.2:4 6:2:30];
% c_prop = [0:.2:4 6:2:30];
% 
% % -- Fixed NM MODEL PARAMETERS --
% mu = 1.3e-7;
% kb = 2.6e-7;
% t_f=mu/kb;
% t_n = 0.01; t_m = .1; %timescales for length, neural, and muscule activity
% a = 1; I = 0; %neural voltage model param, AVB input bias current
% sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
% sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
% dt=1e-3; % time step size
% 
% period = zeros(size(c_ma,2), size(c_prop,2));
% curv_amp = zeros(size(c_ma,2), size(c_prop,2));
% for jj = 1:size(c_ma,2)
%     for kk = 1:size(c_prop,2)
%         
%     % ----  I. FIND PERIODIC ORBIT  ----
%     [ X, period(jj,kk) ] = single_oscillator_LC( dt, c_ma(jj), c_prop(kk), t_f, t_n, t_m, a, I, sigma );
%     max_amp = max(X(:,1))
%     curv_amp(jj,kk) = max_amp;
%     end
% end

% save('oscillator_properties_cees_loop.mat')
load('oscillator_properties_cees_loop.mat');

fig1 = figure(1); clf;
surfc(c_ma, c_prop, period'); view(2);
% title('Period');
xlabel('c_m'); ylabel('c_p'); 
hcb=colorbar();
set(get(hcb,'label'),'string','Period')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig1,'figs/osc_prop_figs/cm_cp_period.png')

fig2 = figure(2); clf;
surfc(c_ma, c_prop, curv_amp'); view(2);
% title('Max Curvature');
xlabel('c_m'); ylabel('c_p');
hcb=colorbar();
set(get(hcb,'label'),'string','Max Curvature')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig2,'figs/osc_prop_figs/cm_cp_curv.png')

fig3 = figure(3); clf;
surfc(c_ma, c_prop, log10(period)'); view(2);
% title('log_{10} Period');
xlabel('c_m'); ylabel('c_p'); 
hcb=colorbar();
set(get(hcb,'label'),'string','log_{10} Period')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig3,'figs/osc_prop_figs/cm_cp_log10period.png')

fig4 = figure(4); clf;
surfc(c_ma, c_prop, log10(curv_amp)'); view(2);
% title('log_{10} Max Curvature');
xlabel('c_m'); ylabel('c_p');
hcb=colorbar();
set(get(hcb,'label'),'string','log_{10} Max Curvature')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig4,'figs/osc_prop_figs/cm_cp_log10curv.png')



%loop over tm, tb
t_m = logspace(-1.3, -.6021,20);
t_b = logspace(-2.301, 0, 20);

% -- Fixed NM MODEL PARAMETERS --
t_n = 0.01; %timescales for neural activity
a = 1; I = 0; %neural voltage model param, AVB input bias current
sigma = @(x) (1/2)*tanh(x-2); %thresholding function for muscle activity into tension
sigma_prime = @(x) (1/2)*(sech(x-2)).^2; %derivative of sigma(x)
dt=1e-3; % time step size
c_ma = 10; c_prop = 1;

period = zeros(size(t_m,2), size(t_b,2));
curv_amp = zeros(size(t_m,2), size(t_b,2));
for jj = 1:size(t_b,2)
    for kk = 1:size(t_m,2)
        
    % ----  I. FIND PERIODIC ORBIT  ----
    [ X, period(jj,kk) ] = single_oscillator_LC( dt, c_ma, c_prop, t_b(jj), t_n, t_m(kk), a, I, sigma );
    max_amp = max(X(:,1));
    curv_amp(jj,kk) = max_amp;
    end
end

fig5 = figure(5); clf;
surfc(t_m, t_b, period'); view(2);
% title('Period');
xlabel('\tau_m'); ylabel('\tau_b');
hcb=colorbar();
set(get(hcb,'label'),'string','Period')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig5,'figs/osc_prop_figs/tm_tb_period.png')


fig6 = figure(6); clf;
surfc(t_m, t_b, curv_amp'); view(2);
% title('Max Curvature');
xlabel('\tau_m'); ylabel('\tau_b'); 
hcb=colorbar();
set(get(hcb,'label'),'string','Max Curvature')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.05;
set(gca,'OuterPosition',a);
saveas(fig6,'figs/osc_prop_figs/tm_tb_curv.png')

fig7 = figure(7); clf;
surfc(t_m, t_b, log10(period)'); view(2);
% title('log_{10} Period');
xlabel('\tau_m'); ylabel('\tau_b');
hcb=colorbar();
set(get(hcb,'label'),'string','log_{10} Period')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)+0.01;
set(gca,'OuterPosition',a);
saveas(fig7,'figs/osc_prop_figs/tm_tb_log10period.png')

fig8 = figure(8); clf;
surfc(t_m, t_b, log10(curv_amp)'); view(2);
% title('log_{10} Max Curvature');
xlabel('\tau_m'); ylabel('\tau_b');
hcb=colorbar();
set(get(hcb,'label'),'string','log_{10} Max Curvature')
set(gca,'FontSize',30);
a = get(gca,'OuterPosition');
a(3) = a(3)-0.01;
set(gca,'OuterPosition',a);
saveas(fig8,'figs/osc_prop_figs/tm_tb_log10curv.png')

