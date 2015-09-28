clear all;
close all;
clc;

%%
fs = [10,12];
TextSizes.DefaultAxesFontSize = fs(2);
TextSizes.DefaultTextFontSize = fs(2);
set(0,TextSizes);

%% Options
plot_options.plot_local = 'false';
plot_options.PL_color = 0.2*[1,1,1];
plot_options.opt_color = 0.2*[1,1,1];
plot_options.threshold_color = 0*[1,1,1];
plot_options.boundary.mark = 0;
plot_options.labels.y_always = true;
plot_options.labels.y_name = 'posterior';
plot_options.legend.box = 'on';
plot_options.legend.position = [0.1,0.1,0.2,0.2];
plot_options.MS.only_optimum = true;
plot_options.MS.name_conv = 'optimum';
plot_options.S.name = 'samples';
plot_options.P.name = 'profiles';
plot_options.A.plot_type = 0;
plot_options.CL.plot_type = 1;

plot_options_2D.plot_local = 'false';
plot_options_2D.PL_color = 0.2*[1,1,1];
plot_options_2D.opt_color = 0.2*[1,1,1];
plot_options_2D.threshold_color = 0*[1,1,1];
plot_options_2D.boundary.mark = 0;
plot_options_2D.labels.y_always = true;
plot_options_2D.labels.y_name = 'posterior';
plot_options_2D.legend.box = 'on';
plot_options_2D.legend.position = [0.02,0.9,0.07,0.07];
plot_options_2D.MS.only_optimum = true;
plot_options_2D.MS.name_conv = 'optimum';
plot_options_2D.S.name = 'samples';
plot_options_2D.P.name = 'profiles';
plot_options_2D.A.plot_type = 0;
plot_options_2D.CL.plot_type = 1;
plot_options_2D.fontsize.tick = 10;


%% M1
% Load data
load('./project/results/model_1__alpha__beta__a0_delta__grad_on__PL.mat'); parameters_P = parameters;
load('./project/results/model_1__alpha__beta__a0_delta__grad_on__MCMC.mat');
parameters.P = parameters_P.P;

% Reordering
ind = [1:8,11,12,9,10];
parameters_new = parameters;
parameters_new.MS.par = parameters.MS.par(ind,:);
parameters_new.MS.par(8,:) = 1 - parameters_new.MS.par(8,:);
parameters_new.S.par = parameters.S.par(ind,:);
parameters_new.S.par(8,:) = 1 - parameters_new.S.par(8,:);
for i = 1:parameters_new.number
    parameters_new.P(i).par = parameters.P(ind(i)).par(ind,:);
    parameters_new.P(i).logPost = parameters.P(ind(i)).logPost;
    parameters_new.P(i).R = parameters.P(ind(i)).R;
    parameters_new.P(i).par(8,:) = 1 - parameters_new.P(i).par(8,:);
end
parameters = parameters_new;

% Options
plot_options.subplot_size_1D = [4,4];
plot_options.subplot_indexing_1D = [1,3,5:8,10,14,11,12,15,16];

% Plot
plotParameterUncertainty(parameters,'1D',[],[],plot_options);

% Refinement
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),1); set(gca,'xtick',[-0.233,-0.23,-0.227]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),3); set(gca,'xtick',[-1.14,-1.10,-1.06,-1.02]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),8); set(gca,'xtick',[-0.45,-0.44,-0.43]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),11); xlim([6.35,7]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),15); xlim([6.35,7]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),12); xlim([-0.85,-0.3]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),16); xlim([-0.85,-0.3]);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 39 40*0.42])
print('-depsc2','figures/figS7');

%% Convergence
figure;
plot(1:length(parameters.S.logPost),parameters.S.logPost,'k.');
xlabel('sample');
ylabel('log-posterior');

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 30 10])
print('-depsc2','figures/figS6a');

%% 2D
% Plot
parameters.S.par = parameters.S.par(:,1:round(size(parameters.S.par,2)/1000):end);
plotParameterUncertainty(parameters,'2D',[],[],plot_options_2D);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 50 50])
print('-depsc2','figures/figS10');

%% M2
% Load data
load('./project/results/model_2__alpha_a__beta__a0_delta__grad_on__PL.mat'); parameters_P = parameters;
load('./project/results/model_2__alpha_a__beta__a0_delta__grad_on__MCMC.mat');
parameters.P = parameters_P.P;

% Reordering
ind = [1:10,13,14,11,12];
parameters_new = parameters;
parameters_new.MS.par = parameters.MS.par(ind,:);
parameters_new.MS.par(10,:) = 1 - parameters_new.MS.par(10,:);
parameters_new.S.par = parameters.S.par(ind,:);
parameters_new.S.par(10,:) = 1 - parameters_new.S.par(10,:);
for i = 1:parameters_new.number
    parameters_new.P(i).par = parameters.P(ind(i)).par(ind,:);
    parameters_new.P(i).logPost = parameters.P(ind(i)).logPost;
    parameters_new.P(i).R = parameters.P(ind(i)).R;
    parameters_new.P(i).par(10,:) = 1 - parameters_new.P(i).par(10,:);
end
parameters = parameters_new;

% Options
plot_options.subplot_size_1D = [4,4];
plot_options.subplot_indexing_1D = [1:8,10,14,11,12,15,16];

% Plot
plotParameterUncertainty(parameters,'1D',[],[],plot_options);

% Refinement
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),11); xlim([6.35,7.4]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),15); xlim([6.35,7.4]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),12); xlim([-0.8,-0.1]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),16); xlim([-0.8,-0.1]);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 39 40*0.42])
print('-depsc2','figures/figS9');

%% Convergence
figure;
plot(1:length(parameters.S.logPost),parameters.S.logPost,'k.');
xlabel('sample');
ylabel('log-posterior');

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 30 10])
print('-depsc2','figures/figS6b');

%% 2D
% Plot
parameters.S.par = parameters.S.par(:,1:round(size(parameters.S.par,2)/1000):end);
plotParameterUncertainty(parameters,'2D',[],[],plot_options_2D);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 50 50])
print('-depsc2','figures/figS11');

%% M16
% Load data
load('./project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on__PL.mat'); parameters_P = parameters;
load('./project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on__MCMC__SampleSize_500000.mat');
parameters.P = parameters_P.P;

% Truncation of warm-up phase
parameters.S.par = parameters.S.par(:,150001:end);
parameters.S.logPost = parameters.S.logPost(150001:end);

% Options
plot_options.subplot_size_1D = [9,4];
plot_options.subplot_indexing_1D = [1:2,5:11,13:14,17:24,25:28,30,34,31,32,35,36];

% Plot
plotParameterUncertainty(parameters,'1D',[],[],plot_options);

% Refinement
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),1); set(gca,'xtick',[-0.318,-0.314,-0.310,-0.306]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),31); xlim([6.35,6.85]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),35); xlim([6.35,6.85]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),32); xlim([-0.85,-0.3]);
subplot(plot_options.subplot_size_1D(1),plot_options.subplot_size_1D(2),36); xlim([-0.85,-0.3]);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 39 40])
print('-depsc2','figures/figS9');

%% Convergence
figure;
plot(1:length(parameters.S.logPost(1:10:end)),parameters.S.logPost(1:10:end),'k.');
xlabel('sample');
ylabel('log-posterior');

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 30 10])
print('-depsc2','figures/figS6c');

%% 2D
fs = [6,8];
TextSizes.DefaultAxesFontSize = fs(2);
TextSizes.DefaultTextFontSize = fs(2);
set(0,TextSizes);

plot_options_2D.fontsize.tick = 6;

% Plot
parameters.S.par = parameters.S.par(:,1:round(size(parameters.S.par,2)/1000):end);
plotParameterUncertainty(parameters,'2D',[],[],plot_options_2D);

% Print
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 50 50])
print('-depsc2','figures/figS12');

