clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

%% PROCESS
% X_1 -> X_2, rate = k_1*[X_1]
% X_2 -> X_1, rate = k_2*[X_2]
% Y = X_2

%% DATA
t = [0:0.2:10]';
ym = [   0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
    1.8309
    3.3559
    4.6091
    5.4235
    5.9757
    6.6298
    7.0080
    7.8280
    7.5852
    7.9247
    7.8363
    8.0107
    7.7077
    7.5316
    7.4208
    7.5734
    7.3197
    7.1489
    7.1987
    6.8493
    6.6425
    6.6268
    6.1223
    6.1078
    5.9242
    5.6034
    5.4618
    5.1281
    4.9489
    4.8930
    4.7747
    4.7750
    4.3095
    4.2211
    4.0416
    3.7485
    3.7164
    3.4799
    3.5286
    3.2785];

%% DEFINITION OF PARAMETER ESTIMATION PROBLEM
% Parameters
parameters.name = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
parameters.min = [           -2;-5;-5;-5;-2];
parameters.max = [log10(max(t)); 5; 5; 5; 2];
parameters.number = 5;

% Properties
properties.name = {'log_{10}(t_0)';'log_{10}(k_{TL}*m_0)';'log_{10}(\beta)';'log_{10}(\delta)';'log_{10}(\sigma)'};
properties.function = {@(theta) propertyFunction_theta(theta,1),...
                       @(theta) propertyFunction_theta(theta,2),...
                       @(theta) propertyFunction_theta(theta,3),...
                       @(theta) propertyFunction_theta(theta,4),...
                       @(theta) propertyFunction_theta(theta,5)};
properties.min = [           -2;-5;-5;-5;-2];
properties.max = [log10(max(t)); 5; 5; 5; 2];
properties.number = length(properties.name);

% Log-posterior function
options.obj_type = 'log-posterior';
logP = @(theta) logP__T(theta,t,ym);

% Log-likelihood and log-prior sperately for multi-chain sampling
logL_and_logPrior = @(theta) logL_and_logPrior__T(theta,t,ym);

%% MULTI-START LOCAL OPTIMIZATION
% Options
options.n_starts = 30;
options.comp_type = 'sequential'; options.mode = 'visual';
% options.comp_type = 'parallel'; options.mode = 'silent'; % n_workers = 10;
% options.save = 'true'; options.foldername = 'results';

% Open matlabpool
if strcmp(options.comp_type,'parallel')
    matlabpool(n_workers);
end

% Optimization
parameters = getMultiStarts(parameters,logP,options);

%% VISUALIZATION: FIT
if strcmp(options.mode,'visual')
    % Simulation
    tsim = linspace(t(1),t(end),100);
    ysim = sim__T(10.^parameters.MS.par(:,1),tsim);

    % Plot: Fit
    figure;
    plot(t,ym,'bo'); hold on;
    plot(tsim,ysim,'r-'); 
    xlabel('time t');
    ylabel('output y');
    legend('data','fit');
end

%% Profile likelihood calculation -- Parameters
parameters = getParameterProfiles(parameters,logP,options);

%% Single-chain Markov chain Monte-Carlo sampling -- Parameters
% options.sampling_scheme = 'DRAM';
% options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'MALA'; options.w_hist = 0;
% options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'MALA'; options.w_hist = 0.5;
% options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'MALA'; options.w_hist = 1;
options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'difference'; options.AM.memory_length = 10*parameters.number;
% options.sampling_scheme = 'single-chain'; options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'position'; options.AM.memory_length = inf;

options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;

parameters = getParameterSamples(parameters,logP,options);

%% Multi-chain Markov chain Monte-Carlo sampling -- Parameters
options.sampling_scheme = 'multi-chain';
options.MC.n_temps = 20;
options.MC.exp_temps = 5;

%options.proposal_scheme = 'MALA'; options.w_hist = 0;
%options.proposal_scheme = 'MALA'; options.w_hist = 0.5;
%options.proposal_scheme = 'AM'; options.AM.adaption_scheme = 'position'; options.AM.memory_length = inf;
options.proposal_scheme = 'AM';
options.AM.adaption_scheme = 'difference';
options.AM.memory_length = 2e2;
options.MC.swapStategy = 'PTEE';

options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;

parameters = getParameterSamples(parameters,logL_and_logPrior,options);

% Visualization of sampling results
options_plot.S.plot_type = 1;
options_plot.S.PT.plot_type = 0; 
options_plot.S.PT.ind = 1:2:options.MC.n_temps;
plotParameterUncertainty(parameters,[],[],[],options_plot);

%% Confidence interval evaluation -- Parameters
alpha = [0.9,0.95,0.99];
parameters = getParameterConfidenceIntervals(parameters,alpha);

%% Evaluation of properties for multi-start local optimization results -- Properties
properties = getPropertyMultiStarts(properties,parameters,options);

%% Profile likelihood calculation -- Properties
properties = getPropertyProfiles(properties,parameters,logP,options);

%% Evaluation of properties for sampling results -- Properties
properties = getPropertySamples(properties,parameters,options);

%% Confidence interval evaluation -- Properties
properties = getPropertyConfidenceIntervals(properties,alpha);

%% Comparison of calculated parameter profiles
if strcmp(options.mode,'visual')
    % Open figure
    figure
    
    % Loop: parameters
    for i = 1:parameters.number
        subplot(ceil(parameters.number/ceil(sqrt(parameters.number))),ceil(sqrt(parameters.number)),i);
        plot(parameters.P(i).par(i,:),parameters.P(i).R,'bx-'); hold on;
        plot(properties.P(i).prop,properties.P(i).R,'r-o');
        xlabel(properties.name{i});
        ylabel('likelihood ratio');
        if i == 1
            legend('unconst. opt. (= standard)','unconst. op. (= new)');
        end
    end
end

%
if strcmp(options.comp_type,'parallel')
    matlabpool('close');
end
