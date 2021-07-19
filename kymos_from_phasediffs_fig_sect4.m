 %--------- show kymographs of traveling wave at 3 fluid viscosities
addpath('./src');
little_gamma = 1.0e+04.*[ 0.0001    0.0010    0.0100    0.0215    0.0464    0.1000    1.0000    2.8003]';

load('sixbox_osc_vs_full_gamma_loop_data.mat');

phase_model_phasediffs = mod(phase_model_eq2,1);
phasediffs1 = phase_model_phasediffs(1,:);
phasediffs2 = phase_model_phasediffs(9,:);
phasediffs3 = phase_model_phasediffs(end,:);

fig2 = figure(2); clf;
subplot(1,3,1); plot(1:5,phasediffs1,'o','MarkerFaceColor','b','MarkerSize',10); 
ylim([0.5 1]); ylabel('\phi_i');
xlabel('i'); set(gca,'FontSize',20,...
    'XTick',[1 2 3 4 5]);
xlim([0 6]);
subplot(1,3,2); plot(1:5,phasediffs2,'o','MarkerFaceColor','b','MarkerSize',10); 
ylim([0.5 1]); ylabel('\phi_i');
xlabel('i'); set(gca,'FontSize',20,...
    'XTick',[1 2 3 4 5]); xlim([0 6]);
subplot(1,3,3); plot(1:5,phasediffs3,'o','MarkerFaceColor','b','MarkerSize',10); 
ylim([0.5 1]); ylabel('\phi_i');
xlabel('i'); set(gca,'FontSize',20,...
    'XTick',[1 2 3 4 5]); xlim([0 6]);

fig1 = figure(1); clf;
dt0 = 0.1*t_n;
%traveling wave (low \gamma)
[X_LC,period] = single_oscillator_LC( dt0, c_ma, c_prop, t_f, t_n, t_m, a, I, sigma );
Xt = phasediffs_to_full_timetrace( dim, phasediffs1, X_LC, period, dt0, gridsz );
subplot(1,3,1); 
Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([round(size(Kappa,2)/2) size(Kappa,2)]);
colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
ylabel('time \rightarrow'); 
xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20,...
    'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
title('\mu_f = 1 mPa s'); %set(gca,'FontSize',20)

%traveling wave (med \gamma)
Xt = phasediffs_to_full_timetrace( dim,phasediffs2, X_LC, period, dt0, gridsz );
subplot(1,3,2); 
Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([round(size(Kappa,2)/2) size(Kappa,2)]);
colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
ylabel('time \rightarrow'); 
xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20,...
    'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
title('\mu_f = 348 mPa s'); %set(gca,'FontSize',20)

%traveling wave (high \gamma)
Xt = phasediffs_to_full_timetrace( dim, phasediffs3, X_LC, period, dt0, gridsz );
subplot(1,3,3); 
Kappa = Xt(1:gridsz*dim,:);
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlim([1 gridsz*dim+1]); ylim([round(size(Kappa,2)/2) size(Kappa,2)]);
colormap(blueblackred); %caxis([2*min(Kappa(1,:)) 2*max(Kappa(1,:))]); colorbar();
ylabel('time \rightarrow'); 
xlabel('body coordinate'); set(gca,'YTickLabel',[],'FontSize',20,...
    'XTick',[1 7], 'XTickLabel', ['0'; 'L';]);
title('\mu_f = 2.8 \times 10^4 mPa s'); %set(gca,'FontSize',20)
    

