% Estimation for M2 for increased upper bound for 

clear all;
close all;
clc;

%% DEFINITION OF MODEL
model_2__alpha_a__beta__a0_delta;

%% LOAD ESTIMATION RESULTS
load('./project/results/model_2__alpha_a__beta__a0_delta__grad_on',...
    'D','a_sim','t_sim','options_logL','options','parameters')

%% MODIFY BOUNDS
parameters_wb.sym = parameters.sym;
parameters_wb.name = parameters.name;
parameters_wb.number = parameters.number;
parameters_wb.min = parameters.min;
parameters_wb.max = parameters.max; parameters_wb.max(2) = 4;
parameters_wb.constraints = parameters.constraints;
parameters_wb.guess = parameters.MS.par(:,1);

%% OPTIMIZATION
% Options
options.n_starts = 1;
options.comp_type = 'sequential';
options.mode = 'visual';
options.fh = figure;
options_logL.fh = figure;

% Timer for visulaization of fit
global tp
tp = clock;

% Optimization
parameters_wb = getMultiStarts(parameters_wb,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options);

%% VISUALIZATION OF OLD AND NEW RESULTS
% n_alpha < 10^2
[Sim] = CPsimulateDALSP(M,parameters.MS.par(:,1),D.t,a_sim,...
                        [D(1).bins(1,1); D(1).bins(:,2)]',...
                        options_logL.simulation);
plotProliferationAssay(Sim,D);

% print
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 26 14]);
print('-depsc2','-r1200',['./figures/RevisionFig_WiderBounds_a_normal']);


% n_alpha < 10^4
[Sim_wb] = CPsimulateDALSP(M,parameters_wb.MS.par(:,1),D.t,a_sim,...
                        [D(1).bins(1,1); D(1).bins(:,2)]',...
                        options_logL.simulation);
plotProliferationAssay(Sim_wb,D);

% print
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 26 14]);
print('-depsc2','-r1200',['./figures/RevisionFig_WiderBounds_b_wider']);

%% VISUALIZATION WITH STEP FUNCTION
a = 0:0.001:1.2;

figure;

% n_alpha < 10^2
n_alpha = 10^parameters.MS.par(2,1);
K_alpha = 10^parameters.MS.par(1,1);
plot(a,(a.^n_alpha)./(K_alpha^n_alpha + a.^n_alpha),'-','color',[0,0,1]); hold on;
plot(a,a > K_alpha,'-','color',[1,0,0]);

xlabel('cell age [d]');
ylabel('function value');

legend('estimated age-dependence','Heaviside step function')

% print
set(gcf,'PaperUnits','centimeters', 'PaperPosition',[0 0 14 7]);
print('-depsc2','-r1200',['./figures/RevisionFig_stepfunction']);

