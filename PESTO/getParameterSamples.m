% getParameterSamples.m performs adaptive MCMC sampling of the posterior
%   distribution using the DRAM tooparameters.minox. The main porpuse of 
%   this routine is to provide a nice interface.  
%
% USAGE:
% ======
% [...] = getParameterSamples(parameters,logPosterior)
% [...] = getParameterSamples(parameters,logPosterior,options)
% [parameters] = getParameterSamples(...)
% [parameters,fh_logPost_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace] = getParameterSamples(...)
% [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis] = getParameterSamples(...)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least:
%   .number ... number of parameter
%   .ml .. maximum likelihood estimate
%   .min ... lower bound for parameter values       
%   .max ... upper bound for parameter values       
% logPosterior ... log-posterior of model as function of the parameters.
% options ... options of algorithm
%   .nsimu_warmup ... length of MCMC warm-up run (default = 1e4).
%   .nsimu_run ... length of MCMC run(default = 5e4).
%   .algorithm ... MCMC sampling scheme (default = 'dram')
%   .qcov ... initial covariance matrix for MCMC sampling
%       (default = 0.001*eye(parameters.number)).
%   .adaptint ... number of function evaluations between adaptation of
%       the MCMC transition kernel (default = 20*parameters.number).
%   .plot ... visualization of the results after the computation (default = 'true').
%   .plot_options ... plot options:
%       .interval ... method uses to determine plot bounds (default = 'dynamic').
%       .hold_on ... conserve of current plot content (default = 'false').
%   .fh_logPost_trace ... figure handle for log-posterior trace plot.
%   .fh_par_trace ... figure handle for parameter trace plots.
%   .fh_par_dis ... figure handle for the parameter distribution plot.
%   .rng ... initialization of random number generator (default = 0).
%       = any ral number r => random generator is initialized with r.
%       = [] ... random number generator is not initialized.
%       (Initializing the random number generator can be helpfult to
%       understand problems.)
%   .plot_options ... plot options for plotPropertyProfiles.m.
%
% Outputs:
% ========
% parameters ... updated parameter object containing:
%   .S ... parameter and posterior sample.
%       .logPost ... log-posterior function along chain
%       .par  ... parameters along chain
% fh_logPost_trace .. figure handle for log-posterior trace
% fh_par_trace .. figure handle for parameter traces
% fh_par_dis .. figure handle for parameter distribution
%
% 2012/07/11 Jan Hasenauer
% 2015/04/29 Jan Hasenauer

% function [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis] = getParameterSamples(parameters,logLikelihood,options)
function [parameters,fh_logPost_trace,fh_par_trace,fh_par_dis_1D,fh_par_dis_2D] = getParameterSamples(varargin)


%% Check and assign inputs
if nargin >= 2
    parameters = varargin{1};
    objective_function = varargin{2};
else
    error('getParameterSamples requires at least three inputs.')
end

% CHANGE JH:
% Check dimension of parameters.min and parameters.max
parameters.min = parameters.min(:);
parameters.max = parameters.max(:);
if ~isfield(parameters,'number')
    parameters.number = length(parameters.min);
end
if    (length(parameters.min) ~= parameters.number) ...
   || (length(parameters.max) ~= parameters.number)
    error('Dimension of parameters.min, parameters.max and parameters.number is inconsistent.');
end

% Set options
options.comp_type = 'sequential'; % 'parallel';
options.obj_type = 'log-posterior'; % 'negative log-posterior'
options.mode = 'visual'; % 'text','silent'
options.save = 'false'; % 'true'
options.rng = 0;
options.thinning = 10;
options.nsimu_warmup = 1e3;
options.nsimu_run    = 1e4;
options.sampling_scheme = 'DRAM';
options.DRAM.algorithm = 'dram';

options.proposal_scheme = 'MALA'; % 'AM'
options.theta_0 = zeros(parameters.number,1);
options.Sigma_0 = 0.001*eye(parameters.number);

% MALA options
options.MALA.min_regularisation = 1e-6; % minimal regularistion for hessian matrix
options.MALA.w_hist = 0; % interpolation between MALA and AM proposal

% AM options
options.AM.min_regularisation = 1e-6; % minimal regularistion for covariance matrix
options.AM.init_memory_length = 20*parameters.number;
options.AM.adaption_interval = 10*parameters.number;
options.AM.min_acc = 0.15;
options.AM.max_acc = 0.30;
options.AM.adap_sigma_scale = 0.8;
% Finite memory:
options.AM.adaption_scheme = 'difference';
options.AM.memory_length = 10*parameters.number;

options.AM.adapt_temperatures = 1; % -> move to MC
options.AM.lacki15_tune_alpha = 0.9; % unify naming
options.AM.start_iter_temp_adaption = 20*parameters.number;
options.AM.proposal_scaling_scheme = 'Lacki15';

% MC options
options.MC.swapStrategy = 'PTEE';
options.MC.n_temps = 10;


options.SCMC.n_proposals = 10;

% % Classical:
% options.AM.adaption_scheme = 'position';
% options.AM.memory_length = inf;

options.plot_options.interval = 'dynamic';
options.fh.logPost_trace = [];
options.fh.par_trace = [];
options.fh.par_dis_1D = [];
options.fh.par_dis_2D = [];

if isfield(parameters,'MS')
    options.theta_0 = parameters.MS.par(:,1);
    if isfield(parameters.MS,'hessian')
        options.Sigma_0 = inv(parameters.MS.hessian(:,:,1));
    end
else
    if isfield(parameters,'guess')
        options.theta_0 = parameters.guess;
    end
end

if nargin == 3
    options = setdefault(varargin{3},options);
end

%% Initialization of random number generator
if ~isempty(options.rng)
    rng(options.rng);
end

%% Initialization and figure generation
fh_logPost_trace = [];
fh_par_trace = [];
fh_par_dis_1D = [];
switch options.mode
    case 'visual'
        % logL trace
        if isempty(options.fh.logPost_trace)
            fh_logPost_trace = figure;
        else
            fh_logPost_trace = figure(options.fh.logPost_trace);
        end
        % parameter traces
        if isempty(options.fh.par_trace)
            fh_par_trace = figure;
        else
            fh_par_trace = figure(options.fh.par_trace);
        end
        % parameter distribution
        if isempty(options.fh.par_dis_1D)
            fh_par_dis_1D = figure;
        else
            fh_par_dis_1D = figure(options.fh.par_dis_1D);
        end
        if isempty(options.fh.par_dis_2D)
            fh_par_dis_2D = figure;
        else
            fh_par_dis_2D = figure(options.fh.par_dis_2D);
        end
    case 'text'
        fprintf(' \nSampling:\n=========\n');
    case 'silent' % no output
end



%% Selection of sampling proceedure
switch options.sampling_scheme
    
    %% DRAM
    case 'DRAM'
        % This section provides the interface to the MATLAB tooparameters.minox for
        % delayed rejection adaptive metropolis sampling developed by 
        % H. Haario et al. (2006), DRAM: Efficient adaptive MCMC,
        % Stat. Comp., 4(16):339-354.
        
        % Model
        for i = 1:parameters.number
            params{i} = {parameters.name{i},options.theta_0(i),parameters.min(i),parameters.max(i)};
        end

        model.ssfun = @(theta,dummi) 2*logPost(theta,objective_function,options.obj_type,'negative');
        model.sigma2 = 1;
        model.N = 1;

        % Options
        options_dram.method      = options.DRAM.algorithm; % adaptation method (mh,am,dr,dram)
        options_dram.qcov        = options.Sigma_0;      % proposal covariance
        options_dram.adaptint    = options.AM.adaption_interval;  % adaptation interval
        options_dram.printint    = 0;  % how often to show info on acceptance ratios
        options_dram.verbosity   = 2;  % how much to show output in Matlab window
        options_dram.updatesigma = 0;  % update error variance
        options_dram.stats       = 0;  % save extra statistics in results
        if strcmp(options.mode,'visual')
            options_dram.waitbar     = 1;  % show garphical waitbar
        else
            options_dram.waitbar     = 0;
        end

        % Warm-up
        options_dram.nsimu = options.nsimu_warmup; % # simulations
        [results] = mcmcrun(model,[],params,options_dram);

        % Sampling
        options_dram.nsimu = options.nsimu_run; % # simulations
        % B: More information
        [results,Theta,~,Obj] = mcmcrun(model,[],params,options_dram,results);

        % Reassignment
        % B: More information!
        parameters.S = results;
        parameters.S.r = -0.5*Obj;
        parameters.S.par = Theta';

        
    case 'single-chain'            

        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.nsimu_run),1);
        j = 0;
        acc = 0;
        
        theta = options.theta_0;
        mu_hist = options.theta_0;
        sigma_hist = options.Sigma_0;
        
        sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive');
                mu = theta;
                Sigma = options.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive');
                if logP < inf
                [mu,Sigma] = getProposal(theta,G,H,options.MALA.min_regularisation,options.MALA.w_hist,...
                    options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end
        
        % Generate Markov chain
        for i = 1:(options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            if mod(i,100) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                             ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
                    case 'text'
                      % B: Added some information and flushing
                      clc
                      disp(str);
                      disp(['Par = ( ' num2str(theta',' %1.2f' ) ' )'])
                    case 'silent' % no output
                end
            end

            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma)';

            % Evaluate objective function
            %B: Transposed to make it more consistent with other notations
            % CHANGE JH: 
            % if (sum(theta_i < parameters.min') + sum(theta_i > parameters.max') == 0)
            if (sum(theta_i < parameters.min) + sum(theta_i > parameters.max) == 0)
                inbounds = 1;
                switch options.proposal_scheme
                    case {'MH','AM'}
                        % Compute log-posterior
                        [logP_i] = logPost(theta_i,objective_function,options.obj_type,'positive');

                        % Update mu and Sigma of proposal
                        mu_i = theta_i;
                        Sigma_i = Sigma;
                    case 'MALA'
                        % Compute log-posterior, gradient and hessian
                        [logP_i,G_i,H_i] = logPost(theta_i,objective_function,options.obj_type,'positive');

                        % Update mu and Sigma of proposal
                        if logP_i < inf
                        [mu_i,Sigma_i] = getProposal(theta_i,G_i,H_i,options.MALA.min_regularisation,options.MALA.w_hist,...
                            mu_hist,sigma_hist,parameters.min,parameters.max);
                        end
                end
            else
                inbounds = 0;
            end

            % Determine acceptance probability
            if (inbounds == 1) && (logP_i < inf)
                % Transition probabilities
                log_p_forward  = logmvnpdf(theta_i,mu  ,Sigma  );
                log_p_backward = logmvnpdf(theta  ,mu_i,Sigma_i);

                % Acceptance probability
                pacc = exp(logP_i - logP + log_p_backward - log_p_forward);
            else
                pacc = 0;
            end

            % Accept or reject
            r = rand;
            if r <= pacc
                acc    = acc + 1;
                theta  = theta_i;
                dtheta = (theta-mu);
                logP   = logP_i;
                mu     = mu_i;
                Sigma  = Sigma_i; % only for MALA relevant
            else
                dtheta = zeros(size(theta));
            end
            
            % Incremental calculation of mean and covariance
            % (with memory length options.AM.memory_length)
            switch options.AM.adaption_scheme
                case 'position'
                    [mu_hist,sigma_hist] = updateStatistics(mu_hist,sigma_hist,theta,max(i,options.AM.init_memory_length),...
                        sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                case 'difference'
                    [sigma_hist] = updateCovariance(sigma_hist,dtheta,max(i,options.AM.init_memory_length),...
                        sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    mu_hist = theta;
            end
            
            % Proposal update
            if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
                if acc/i < options.AM.min_acc
                    sigma_scale = sigma_scale*options.AM.adap_sigma_scale;
                elseif acc/i > options.AM.max_acc
                    sigma_scale = sigma_scale/options.AM.adap_sigma_scale;
                end
                Sigma = sigma_scale*sigma_hist;
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.AM.min_regularisation*eye(parameters.number);
                end
            end
            
            % Store
            if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                j = j + 1;
                parameters.S.par(:,j) = theta;
                parameters.S.logPost(j) = logP;
            end        
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
    case 'single-chain multi-core'            

        n_proposals = options.SCMC.n_proposals;
        useful_total = 0;
        not_useful_total = 0;
        K_set = [];
       
        % Initialization
        parameters.S.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run));
        parameters.S.logPost = nan(length(1:options.thinning:options.nsimu_run),1);
        j = 0;
        acc = 0;
        K = 0;

        theta = options.theta_0;
        mu_hist = options.theta_0;
        sigma_hist = options.Sigma_0;
        
        sigma_scale = 1;
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                [logP] = logPost(theta,objective_function,options.obj_type,'positive');
                mu = theta;
                Sigma = options.Sigma_0;
            case 'MALA'
                [logP,G,H] = logPost(theta,objective_function,options.obj_type,'positive');
                if logP < inf
                [mu,Sigma] = getProposal(theta,G,H,options.MALA.min_regularisation,options.MALA.w_hist,...
                    options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                end
        end
        if isnan(logP) || (logP == -inf)
            error('log-posterior undefined at initial point.');
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end
        
        % Generate Markov chain
        i = 1;
        while i <= (options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            if min(mod([i-K:i],100)) == 0
                str = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                             ' % (acc = ' num2str(100*acc/i,'%.2f') ' % )'];
                switch options.mode
                    case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
                    case 'text', disp(str);
                    case 'silent' % no output
                end
            end

            % Propose new parameter vector
            theta_i = mvnrnd(mu,Sigma,n_proposals)';

            % Evaluate objective function
            for k = 1:n_proposals
                if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max) == 0)
                    inbounds(k) = 1;
                    switch options.proposal_scheme
                        case {'MH','AM'}
                            % Compute log-posterior
                            [logP_i(k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            mu_i(:,k) = theta_i(:,k);
                            Sigma_i(:,:,k) = Sigma;
                        case 'MALA'
                            % Compute log-posterior, gradient and hessian
                            [logP_i(k),G_i,H_i] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            if logP_i(k) < inf
                                [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),G_i,H_i,options.MALA.min_regularisation,options.MALA.w_hist,mu_hist,sigma_hist,parameters.min,parameters.max);
                            end
                    end
                else
                    inbounds(k) = 0;
                end
                
                % Determine acceptance probability
                if (inbounds(k) == 1) && (logP_i(k) < inf)
                    % Transition probabilities
                    log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu       ,Sigma         );
                    log_p_backward(k) = logmvnpdf(theta       ,mu_i(:,k),Sigma_i(:,:,k));

                    % Acceptance probability
                    pacc(k) = exp(logP_i(k) - logP + log_p_backward(k) - log_p_forward(k));
                else
                    pacc(k) = 0;
                end

                % Accept or reject
                if rand <= pacc(k)
                    acc_flag(k) = 1;
                else
                    acc_flag(k) = 0;
                end

            end

            % Assignment
            K = min(find(acc_flag));
            if isempty(K)
                K = n_proposals;
            end
            K_set(end+1) = K;
            useful_total = useful_total + K;
            not_useful_total = not_useful_total + n_proposals - K;
            [K,useful_total,not_useful_total]
            for k = 1:K
                % Assignment of new state
                if acc_flag(k) == 1
                    acc    = acc + 1;
                    theta  = theta_i(:,k);
                    dtheta = (theta-mu);
                    logP   = logP_i(k);
                    mu     = mu_i(:,k);
                    Sigma  = Sigma_i(:,:,k); % only for MALA relevant
                else
                    dtheta = zeros(size(theta));
                end
                            
                % Incremental calculation of mean and covariance
                % (with memory length options.AM.memory_length)
                switch options.AM.adaption_scheme
                    case 'position'
                        [mu_hist,sigma_hist] = updateStatistics(mu_hist,sigma_hist,theta,max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    case 'difference'
                        [sigma_hist] = updateCovariance(sigma_hist,dtheta,max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                        mu_hist = theta;
                    case 'Lacki15'
                        [mu_hist,sigma_hist] = ...
                          updateStatistics_Lacki15(mu_hist,...
                                                   sigma_hist,...
                                                   theta,...
                                                   max(i,options.AM.init_memory_length),...
                                                   sqrt(2)/options.AM.memory_length,...
                                                   options.AM.lacki15_tune_alpha,...
                                                   options.AM.min_regularisation);                        
                end

                % Proposal update
                if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
                    if acc/i < options.AM.min_acc
                        sigma_scale = sigma_scale*options.AM.adap_sigma_scale;
                    elseif acc/i > options.AM.max_acc
                        sigma_scale = sigma_scale/options.AM.adap_sigma_scale;
                    end
                    Sigma = sigma_scale*sigma_hist;
                    switch options.AM.proposal_scaling_scheme
                      case 'original'
                          if acc/i < options.AM.min_acc
                              sigma_scale = sigma_scale*options.AM.adap_sigma_scale;
                          elseif acc/i > options.AM.max_acc
                              sigma_scale = sigma_scale/options.AM.adap_sigma_scale;
                          end                    
                      case 'Lacki15'
                        	sigma_scale = log(sigma_scale)/2 + ...
                            (acc/i - 0.234) / (i+1)^options.AM.lacki15_tune_alpha;
                          sigma_scale = exp(2*sigma_scale);
                    end
                end
                
                % Regularisation
                [~,p] = cholcov(Sigma,0);
                if p ~= 0
                    Sigma = Sigma + options.AM.min_regularisation*eye(parameters.number);
                end
                
                % Store
                if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                    j = j + 1;
                    parameters.S.par(:,j) = theta;
                    parameters.S.logPost(j) = logP;
                end
                
                % Update of counter
                i = i + 1;
            end
            
        end
        % Reduction
        parameters.S.par = parameters.S.par(:,1:j);
        parameters.S.logPost = parameters.S.logPost(1:j);
        
        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end
        
        save('test','K_set');
    case 'multi-chain'

        % Initialization
        beta = linspace(0,1,options.MC.n_temps).^options.MC.exp_temps;
        % B: 0 excluded
        beta = fliplr(linspace(beta(2),1,options.MC.n_temps).^options.MC.exp_temps);
        
        % B:
        outOfBounds = zeros(1,options.MC.n_temps);

        j = 0;
        parameters.S.PT.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
        parameters.S.PT.logPost = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
        
        % B:
        parameters.S.PT.temp_reduction_iters = [];
        parameters.S.PT.temp_reduction_numbers = [];

        acc = zeros(options.MC.n_temps,1);
        % B: allgemeinere Übergänge möglich, daher Matrix
        acc_swap = zeros(options.MC.n_temps,options.MC.n_temps);

        mu = nan(parameters.number,options.MC.n_temps);
        mu_i = nan(parameters.number,options.MC.n_temps);
        mu_hist = repmat(options.theta_0,[1,options.MC.n_temps]);
        Sigma = nan(parameters.number,parameters.number,options.MC.n_temps);
        Sigma_i = nan(parameters.number,parameters.number,options.MC.n_temps);
        
        % B: Init new properties
        parameters.S.PT.acc = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
        parameters.S.PT.acc_swap = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps,options.MC.n_temps);
        parameters.S.PT.sigma_scale = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
        
        
        % B: Possibility to include different starting proposal Sigmas for
        % each chain
        if length(size(options.Sigma_0)) == 3
          sigma_hist = options.Sigma_0;
        else
          sigma_hist = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);
        end

        sigma_scale = ones(options.MC.n_temps,1);
        
        logL = nan(options.MC.n_temps,1);
        logL_i = nan(options.MC.n_temps,1);
        logPrior = nan(options.MC.n_temps,1);
        logPrior_i = nan(options.MC.n_temps,1);
        
        % B: Possibility to include different starting points for each
        % chain
        dummy = size(options.theta_0);
        if min(dummy(1:2)) > 1
          theta = options.theta_0;
        else
          theta = repmat(options.theta_0,[1,options.MC.n_temps]);
        end
        dtheta = zeros(parameters.number,options.MC.n_temps);
        
        % Initialization and testing of starting point
        switch options.proposal_scheme
            case {'MH','AM'}
                % Objective function evaluation
                [logL_0,logPrior_0] = ...
                    logPost_pt(options.theta_0,objective_function,options.obj_type,'positive');

                % Assignment
                logL = repmat(logL_0,[options.MC.n_temps,1]);
                logPrior = repmat(logPrior_0,[options.MC.n_temps,1]);
                mu = repmat(options.theta_0,[1,options.MC.n_temps]);
                Sigma = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);
            case 'MALA'
                for k = 1:options.MC.n_temps
                    % Objective function evaluation
                    [logL(k),logPrior(k),dlogL(:,k),dlogPrior(:,k),ddlogL(:,:,k),ddlogPrior(:,:,k)] = ...
                        logPost_pt(options.theta_0,objective_function,options.obj_type,'positive');
                    
                    % Assignment
                    if logL(k) < inf
                    [mu(:,k),Sigma(:,:,k)] = getProposal(options.theta_0,...
                        beta(k)*dlogL(:,k)+dlogPrior(:,k),beta(k)*ddlogL(:,:,k)+ddlogPrior(:,:,k),...
                        options.MALA.min_regularisation,options.MALA.w_hist,...
                        options.theta_0,options.Sigma_0,parameters.min,parameters.max);
                    end
                end
        end
       
        if isnan(logL(end)) || (logL(end) == -inf)
            error(['log-posterior undefined at initial point or error during '...
                   'objective function evaluation in multi-chain sampler. '   ...
                   '(Probably no second input corresponding to temperature.)']);
        end

        % Initialization of waitbar
        if strcmp(options.mode,'visual')
            h = waitbar(0,['Sampling completed to 0 % (acc = 0 %)']);
        end

        % Generate Markov chain
        istr = 'start';
        for i = 1:(options.nsimu_run+options.nsimu_warmup)
            % Report of progress
            % B:stetigere ausgabe
            if mod(i,1) == 0
                istr = ['Sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
                     ' % (acc = ' num2str(100*acc(1)/i,'%.2f') ' % )'];
%                 switch options.mode
%                     case 'visual', waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,str);
%                     case 'text', disp(istr);
%                     case 'silent' % no output
%                 end
            end
            
            for k = 1:options.MC.n_temps
                % Propose new parameter vector
                theta_i(:,k) = mvnrnd(mu(:,k),Sigma(:,:,k))';

                % Evaluate objetive function
                % B: Transpose inconsitent to other functions -> added ' to
                % JH: removed transpose!!!
                % parameters.min/max
                if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max)) == 0
                    inbounds = 1;
                    switch options.proposal_scheme
                        case {'MH','AM'}
                            % Objective function evaluation
                            % B: Funktionsaufruf Übergabeargumente: 
                            % logPrior war bisher 
                            % eigentlich der Gradient
                            % -> neue Struktur? Inkonsistent zur
                            % Optimierung, besser prior als viertes output
                            % argument bei objective function verwenden (2. 
                            % ist gradient, 3. hessian)
                            [logL_i(k),logPrior_i(k)] = ...
                                logPost_pt(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Update mu and Sigma of proposal
                            mu_i(:,k) = theta_i(:,k);
                            Sigma_i(:,:,k) = Sigma(:,:,k);
                        case 'MALA'
                            % Objective function evaluation
                            [logL_i(k),logPrior_i(k),dlogL_i(:,k),dlogPrior_i(:,k),ddlogL_i(:,:,k),ddlogPrior_i(:,:,k)] = ...
                                logPost_pt(theta_i(:,k),objective_function,options.obj_type,'positive');

                            % Assignment
                            if logL_i(k) < inf
                                logL_i(k)
                                ddlogL_i(:,:,k)
                            [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),...
                                beta(k)*dlogL_i(:,k)+dlogPrior_i(:,k),beta(k)*ddlogL_i(:,:,k)+ddlogPrior_i(:,:,k),...
                                options.MALA.min_regularisation,options.MALA.w_hist,...
                                mu_hist(:,k),sigma_hist(:,:,k),parameters.min,parameters.max);
                            end
                    end
                else
                    inbounds = 0;
                end
                
                % Determine acceptance probability
                if (inbounds == 1) && (logL_i(k) < inf)
                    % Transition probabilities
                    % B: Auch im Fall symmetrischer Proposals sinnvoll?
                    log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu(:,k)  ,Sigma(:,:,k)  );
                    log_p_backward(k) = logmvnpdf(theta(:,k)  ,mu_i(:,k),Sigma_i(:,:,k));

                    % Acceptance probability
                    pacc(k) = exp(  (beta(k)*logL_i(k)+logPrior_i(k)) ...
                                  - (beta(k)*logL(k)  +logPrior(k)  ) ...
                                  + log_p_backward(k) - log_p_forward(k));
                else
                    pacc(k) = 0;
                    outOfBounds(k) = outOfBounds(k) + 1;
                end
                
                % Accept or reject
                if rand <= pacc(k)
                    acc(k)       = acc(k) + 1;
                    theta(:,k)   = theta_i(:,k);
                    dtheta(:,k)  = (theta(:,k)-mu(:,k));
                    logL(k)      = logL_i(k);
                    logPrior(k)  = logPrior_i(k);
                    mu(:,k)      = mu_i(:,k);
                    Sigma(:,:,k) = Sigma_i(:,:,k); % only for MALA relevant
                else
                    dtheta(:,k) = 0;
                end
            end
            
            % Loop: Swaps
            if options.MC.n_temps > 1
              switch options.MC.swapStrategy 

                % B: Alle nächsten Nachbarn werden vorgeschlagen
                case 'all_adjacents'
                  pacc_swap = [];
                  for k = 1:options.MC.n_temps-1

                      pacc_swap = exp((beta(k)-beta(k+1))*(logL(k+1)-logL(k)));

                      if rand <= pacc_swap
                        acc_swap(k,k+1)        = acc_swap(k,k+1) + 1;
                        theta(:,[k,k+1])   = theta(:,[k+1,k]);
                        logL([k,k+1])      = logL([k+1,k]);
                        logPrior([k,k+1])  = logPrior([k+1,k]);
                        mu(:,[k,k+1])      = mu(:,[k+1,k]);
                        switch options.proposal_scheme
                            case 'MALA'
                                Sigma(:,:,[k,k+1]) = Sigma(:,:,[k+1,k]); % only for MALA relevant
                        end
                      end
                  end
                case 'reversed_adjacents'
                  pacc_swap = [];
                  for k = options.MC.n_temps-1:-1:1

                      pacc_swap = exp((beta(k)-beta(k+1))*(logL(k+1)-logL(k)));

                      if rand <= pacc_swap
                        acc_swap(k,k+1)        = acc_swap(k,k+1) + 1;
                        theta(:,[k,k+1])   = theta(:,[k+1,k]);
                        logL([k,k+1])      = logL([k+1,k]);
                        logPrior([k,k+1])  = logPrior([k+1,k]);
                        mu(:,[k,k+1])      = mu(:,[k+1,k]);
                        switch options.proposal_scheme
                            case 'MALA'
                                Sigma(:,:,[k,k+1]) = Sigma(:,:,[k+1,k]); % only for MALA relevant
                        end
                      end
                  end

                case 'all_adjacents_randomized'
                  pacc_swap = [];

                  rn = randperm(options.MC.n_temps-1);

                  for k = 1:options.MC.n_temps-1

                      pacc_swap = exp((beta(rn(k))-beta(rn(k)+1))*...
                                          (logL(rn(k)+1)-logL(rn(k))));

                      if rand <= pacc_swap
                        acc_swap(rn(k),rn(k)+1)        = acc_swap(rn(k),rn(k)+1) + 1;
                        theta(:,[rn(k),rn(k)+1])   = theta(:,[rn(k)+1,rn(k)]);
                        logL([rn(k),rn(k)+1])      = logL([rn(k)+1,rn(k)]);
                        logPrior([rn(k),rn(k)+1])  = logPrior([rn(k)+1,rn(k)]);
                        mu(:,[rn(k),rn(k)+1])      = mu(:,[rn(k)+1,rn(k)]);
                        switch options.proposal_scheme
                            case 'MALA'
                                Sigma(:,:,[rn(k),rn(k)+1]) = Sigma(:,:,[rn(k)+1,rn(k)]); % only for MALA relevant
                        end
                      end
                  end

                case 'all'
                  pacc_swap = [];
                  for k = 1:options.MC.n_temps
                    for l = 1:options.MC.n_temps
                      if k ~= l
                        pacc_swap = exp((beta(k)-beta(l))*(logL(l)-logL(k)));
                      end

                      if rand <= pacc_swap
                        acc_swap(k,l)        = acc_swap(k,l) + 1;
                        theta(:,[k,l])   = theta(:,[l,k]);
                        logL([k,l])      = logL([l,k]);
                        logPrior([k,l])  = logPrior([l,k]);
                        mu(:,[k,l])      = mu(:,[l,k]);
                        switch options.proposal_scheme
                            case 'MALA'
                                Sigma(:,:,[k,l]) = Sigma(:,:,[l,k]); % only for MALA relevant
                        end
                      end
                    end
                  end

                % B: Swap Vorschlag basierend auf Posteriorverhältnis für
                % ALLE Kombinationen möglich  
                case 'PTEE'
                  pacc_swap = [];
                  swap_propose = zeros(options.MC.n_temps,options.MC.n_temps);
                  for k = 1:options.MC.n_temps
                    for m = 1:k-1
                      swap_propose(k,m) = exp(-abs(logL(k)-logL(m)));
                    end
                  end

                  perm = rand(options.MC.n_temps,options.MC.n_temps) < swap_propose;
                  [p1,p2] = find(perm);
                  rn = randperm(length(p1)); 
                  p1 = p1(rn);
                  p2 = p2(rn);

                  for k = 1:length(p1)
                    pacc_swap = exp((beta(p1(k))-beta(p2(k)))...
                                          *(logL(p2(k))-logL(p1(k))));

                    if rand <= pacc_swap
                      acc_swap(p1(k),p2(k)) = acc_swap(p1(k),p2(k)) + 1;
                      theta(:,[p1(k),p2(k)])   = theta(:,[p2(k),p1(k)]);
                      logL([p1(k),p2(k)])      = logL([p2(k),p1(k)]);
                      logPrior([p1(k),p2(k)])  = logPrior([p2(k),p1(k)]);
                      mu(:,[p1(k),p2(k)])      = mu(:,[p2(k),p1(k)]);
                      switch options.proposal_scheme
                          case 'MALA'
                              Sigma(:,:,[p1(k),p2(k)]) = Sigma(:,:,[p2(k),p1(k)]); % only for MALA relevant
                      end
                      dummy1 = p1(k);
                      dummy2 = p2(k);
                      p1(p1 == dummy1) = dummy2;
                      p2(p2 == dummy2) = dummy1;
                    end
                  end
              end
            end
            % Loop: Mean, covariance and proposal update
            for k = 1:options.MC.n_temps
                % Incremental calculation of mean and covariance
                % (with memory length options.AM.memory_length)
                switch options.AM.adaption_scheme
                    case 'position' %B: hier war ein Fehler drin
                        [mu_hist(:,k),sigma_hist(:,:,k)] = updateStatistics(mu_hist(:,k),sigma_hist(:,:,k),theta(:,k),max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                    case 'difference' %B: hier war ein Fehler drin
                        [sigma_hist(:,:,k)] = updateCovariance(sigma_hist(:,:,k),dtheta(:,k),max(i,options.AM.init_memory_length),...
                            sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
                        mu_hist(:,k) = theta(:,k);
                    case 'Lacki15'
                        [mu_hist(:,k),sigma_hist(:,:,k)] = ...
                          updateStatistics_Lacki15(mu_hist(:,k),...
                                                   sigma_hist(:,:,k),...
                                                   theta(:,k),...
                                                   max(i,options.AM.init_memory_length),...
                                                   sqrt(2)/options.AM.memory_length,...
                                                   options.AM.lacki15_tune_alpha,...
                                                   options.AM.min_regularisation);
                end
            
                % Proposal update
                if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
                    switch options.AM.proposal_scaling_scheme
                      case 'original'
                          if acc(k)/i < options.AM.min_acc
                              sigma_scale(k) = sigma_scale(k)*options.AM.adap_sigma_scale;
                          elseif acc(k)/i > options.AM.max_acc
                              sigma_scale(k) = sigma_scale(k)/options.AM.adap_sigma_scale;
                          end
                      case 'Lacki15'
                        	sigma_scale(k) = log(sigma_scale(k))/2 + ...
                            (acc(k)/i - 0.234) / (i+1)^options.AM.lacki15_tune_alpha;
                          sigma_scale(k) = exp(2*sigma_scale(k));
                      case 'none'
                          sigma_scale(k) = 1;
                      otherwise
                        error('You have to specify a proper proposal scaling scheme type!')
                    end
                    
                    Sigma(:,:,k) = sigma_scale(k)*sigma_hist(:,:,k);

                    % Regularisation
                    [~,p] = cholcov(Sigma(:,:,k),0);
                    if p ~= 0
                        Sigma(:,:,k) = Sigma(:,:,k) + options.AM.min_regularisation*eye(parameters.number);
                    end
                end
            end
            
            % B: Adaptive temperature amount reduction
            if options.AM.adapt_temperatures == 1 && options.MC.n_temps > 1
              T = 1./beta;
              T_old = T;
              T(1) = 1;
              xi = zeros(1,options.MC.n_temps-1);
              for l = 1 : options.MC.n_temps - 1
                xi(l) = min(1, exp((beta(l+1)-beta(l))*(logL(l)-logL(l+1))));
                T(l+1) = T(l) + (T(l+1)-T_old(l)) * ...
                          exp((xi(l)-0.234)/(i+1)^options.AM.lacki15_tune_alpha);
              end
              beta = 1./T;
              if i > options.AM.start_iter_temp_adaption
                switch options.AM.proposal_scaling_scheme 
                  case 'Lacki15'
                    dummy = (sigma_scale/2 >= 2.38 / sqrt(length(theta(:,1)))...
                                          * ones(options.MC.n_temps,1));
                  otherwise
                    dummy = (exp(sigma_scale) >= 2.38 / sqrt(length(theta(:,1)))...
                                          * ones(options.MC.n_temps,1));
                end
                temps_old = options.MC.n_temps;
                options.MC.n_temps = min([find(dummy)',options.MC.n_temps]);
                
                % B: reduce everything to smaller size
                if temps_old > options.MC.n_temps
                  
                  Sigma = Sigma(:,:,1:options.MC.n_temps);
                  sigma_hist = sigma_hist(:,:,1:options.MC.n_temps);
                  Sigma_i = Sigma_i(:,:,1:options.MC.n_temps);
                  sigma_scale = sigma_scale(1:options.MC.n_temps);
                  beta = beta(1:options.MC.n_temps);
                  acc = acc(1:options.MC.n_temps);
                  acc_swap = acc_swap(1:options.MC.n_temps,1:options.MC.n_temps);
                  dtheta = dtheta(:,1:options.MC.n_temps);
                  logL = logL(1:options.MC.n_temps);
                  logL_i = logL_i(1:options.MC.n_temps);
                  logPrior = logPrior(1:options.MC.n_temps);
                  logPrior_i = logPrior_i(1:options.MC.n_temps);
                  log_p_backward = log_p_backward(1:options.MC.n_temps);
                  log_p_forward = log_p_forward(1:options.MC.n_temps);
                  mu = mu(:,1:options.MC.n_temps);
                  mu_hist = mu_hist(:,1:options.MC.n_temps);
                  mu_i = mu_i(:,1:options.MC.n_temps);
                  outOfBounds = outOfBounds(1:options.MC.n_temps);
                  pacc = pacc(1:options.MC.n_temps);
                  theta = theta(:,1:options.MC.n_temps);
                  theta_i = theta_i(:,1:options.MC.n_temps);  
                  
%                   lazystr = {'par','logPost','acc','acc_swap','sigma_scale',...
%                              'sigma_hist', 'outOfBoundsCount'};
%                   for b = 1:length(lazystr)
%                       if isfield(parameters.S.PT,lazystr{b})
%                         parameters.S.(['nT_' num2str(options.MC.n_temps)]).(lazystr{b}) = ...
%                             parameters.S.PT.(lazystr{b});
%                       end
%                       parameters.S.PT.(lazystr{b}) = [];
%                   end    
                  
                  parameters.S.PT.temp_reduction_iters(end+1) = i;
                  parameters.S.PT.temp_reduction_numbers(end+1) = options.MC.n_temps;
                  
                end
              end
            end

            
            % Store
            if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
                j = j + 1;
                
                % B: Make sure the resizing works properly - not necessary
                % direktly after warm up
                if ~isempty(parameters.S.PT.par)
                  parameters.S.PT.par = parameters.S.PT.par(:,:,1:options.MC.n_temps);
                  parameters.S.PT.logPost = parameters.S.PT.logPost(:,1:options.MC.n_temps);
                  parameters.S.PT.acc = parameters.S.PT.acc(:,1:options.MC.n_temps);
                  parameters.S.PT.acc_swap = ...
                    parameters.S.PT.acc_swap(:,1:options.MC.n_temps,1:options.MC.n_temps);
                  parameters.S.PT.sigma_scale = ...
                    parameters.S.PT.sigma_scale(:,1:options.MC.n_temps);
                end
                parameters.S.PT.par(:,j,:) = theta;
                parameters.S.PT.logPost(j,:) = logL + logPrior;
                
                % B:more information    
                parameters.S.PT.acc(j,:) = 100*acc(:)/i;
                parameters.S.PT.acc_swap(j,:,:) = 100*(acc_swap(:,:)+acc_swap(:,:)')/i;
                parameters.S.PT.sigma_scale(j,:) = sigma_scale(:);
                parameters.S.PT.sigma_hist = sigma_hist(:,:,:);
                parameters.S.PT.outOfBoundsCount = outOfBounds;
                
                str = num2str(100*acc(1)/i,'%.2f');
                if options.MC.n_temps > 1
                  for k = 2:options.MC.n_temps
                      str = [str ', ' num2str(100*acc(k)/i,'%.2f')];
                  end
                end
                % B:flushed display
                clc
                disp(istr);
                bstr = ['iter: ' num2str(i)...
                    '   #Temp = ' num2str(options.MC.n_temps)...
                    '    lowest inverse Temp = ' num2str(beta(end))];
                disp(bstr);
%                 disp([num2str(i,'%i') ' / ' num2str((options.nsimu_run+options.nsimu_warmup))]);
                disp(['acc = ' str ' %']);
                disp(['PAR: ' num2str(theta(:,1)','%1.2f  ')]);
                
                acc_swap2 = acc_swap + acc_swap';
                str = '';
                for k = 1:options.MC.n_temps-1
                    str = [str num2str(100*(sum(acc_swap2(k,:)))/...
                      i,'%.2f') ', '];
                end
                disp(['acc_swap = ' str ' %']);
            end
        end
        
        % Reduction
        parameters.S.PT.par = parameters.S.PT.par(:,1:j,:);
        parameters.S.PT.logPost = parameters.S.PT.logPost(1:j,:);
        parameters.S.par = parameters.S.PT.par(:,1:j,end);
        parameters.S.logPost = parameters.S.PT.logPost(1:j,end);

        % Close waitbar
        if strcmp(options.mode,'visual')
            close(h)
        end

%         % Check
%         try
%             objective_function(options.theta_0);
%         catch
%             error(['Error during objective function evaluation in multi-chain sampler. '  ...
%                    '(Probably no second input corresponding to temperature.)']);
%         end
%             
%         % Initialization
%         temp = linspace(0,1,options.MC.n_temps).^options.MC.exp_temps;
% 
%         parameters.S.PT.par = nan(parameters.number,length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
%         parameters.S.PT.logPost = nan(length(1:options.thinning:options.nsimu_run),options.MC.n_temps);
%         j = 0;
%         acc = zeros(options.MC.n_temps,1);
%         acc_trans = zeros(options.MC.n_temps-1,1);
% 
%         theta = repmat(options.theta_0,[1,options.MC.n_temps]);
%         mu_hist = repmat(options.theta_0,[1,options.MC.n_temps]);
%         sigma_hist = repmat(options.Sigma_0,[1,1,options.MC.n_temps]);
%         
%         sigma_scale = ones(options.MC.n_temps,1);
%         
%         logP = nan(options.MC.n_temps,1);
%         logP_i = nan(options.MC.n_temps,1);
%         dtheta = nan(parameters.number,options.MC.n_temps);
%         
%         % Initialization and testing of starting point
%         for k = 1:options.MC.n_temps
%             switch options.proposal_scheme
%                 case {'MH','AM'}
%                     [logP(k)] = logPost(theta(:,k),objective_function,options.obj_type,'positive',temp(k));
%                     mu(:,k) = theta(:,k);
%                     Sigma(:,:,k) = options.Sigma_0;
%                 case 'MALA'
%                     [logP(k),G(:,k),H(:,:,k)] = logPost(theta(:,k),objective_function,options.obj_type,'positive',temp(k));
%                     [mu(:,k),Sigma(:,:,k)] = getProposal(theta(:,k),G(:,k),H(:,:,k),options.MALA.min_regularisation,options.MALA.w_hist,...
%                         options.theta_0,options.Sigma_0,parameters.min,parameters.max);
%             end
%         end
%         if isnan(logP(end)) || (logP(end) == -inf)
%             error('log-posterior undefined at initial point.');
%         end
% 
%         % Initialization of waitbar
%         h = waitbar(0,['MCMC sampling completed to 0 % (acc = 0 %)']);
% 
%         % Generate Markov chain
%         for i = 1:(options.nsimu_run+options.nsimu_warmup)
%             % Report of progress
%             if mod(i,100) == 0
%                 waitbar(i/(options.nsimu_run + options.nsimu_warmup),h,...
%                     ['MCMC sampling completed to ' num2str(100*i/(options.nsimu_run + options.nsimu_warmup),'%.2f')...
%                      ' % (acc = ' num2str(100*acc(end)/i,'%.2f') ' % )']);
%             end
% 
%             for k = 1:options.MC.n_temps
%                 % Propose new parameter vector
%                 theta_i(:,k) = mvnrnd(mu(:,k),Sigma(:,:,k))';
% 
%                 % Check bounds
%                 if (sum(theta_i(:,k) < parameters.min) + sum(theta_i(:,k) > parameters.max) == 0)
%                     switch options.proposal_scheme
%                         case {'MH','AM'}
%                             % Compute log-posterior
%                             [logP_i(k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',temp(k));
% 
%                             % Update mu and Sigma of proposal
%                             mu_i(:,k) = theta_i(:,k);
%                             Sigma_i(:,:,k) = Sigma(:,:,k);
%                         case 'MALA'
%                             % Compute log-posterior, gradient and hessian
%                             [logP_i(k),G_i(:,k),H_i(:,:,k)] = logPost(theta_i(:,k),objective_function,options.obj_type,'positive',temp(k));
% 
%                             % Update mu and Sigma of proposal
%                             [mu_i(:,k),Sigma_i(:,:,k)] = getProposal(theta_i(:,k),G_i(:,k),H_i(:,:,k),options.MALA.min_regularisation,options.MALA.w_hist,...
%                                 mu_hist(:,k),sigma_hist(:,:,k),parameters.min,parameters.max);
%                     end
% 
%                     % Transition probabilities
%                     log_p_forward(k)  = logmvnpdf(theta_i(:,k),mu(:,k)  ,Sigma(:,:,k)  );
%                     log_p_backward(k) = logmvnpdf(theta(:,k)  ,mu_i(:,k),Sigma_i(:,:,k));
% 
%                     % Acceptance probability
%                     pacc(k) = exp(logP_i(k) - logP(k) + log_p_backward(k) - log_p_forward(k));
%                 else
%                     pacc(k) = 0;
%                 end
% 
%                 % Accept or reject
%                 r(k) = rand;
%                 if r(k) <= pacc(k)
%                     acc(k)       = acc(k) + 1;
%                     theta(:,k)   = theta_i(:,k);
%                     dtheta(:,k)  = (theta(:,k)-mu(:,k));
%                     logP(k)      = logP_i(k);
%                     mu(:,k)      = mu_i(:,k);
%                     Sigma(:,:,k) = Sigma_i(:,:,k); % only for MALA relevant
%                 else
%                     dtheta(:,k) = 0;
%                 end
%                 
%                 % Incremental calculation of mean and covariance
%                 % (with memory length options.AM.memory_length)
%                 switch options.AM.adaption_scheme
%                     case 'position'
%                         [mu_hist(:,k),sigma_hist(:,:,k)] = updateStatistics(mu_hist(:,k),sigma_hist(:,:,k),theta(:,k),max(i,options.AM.init_memory_length),...
%                             sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
%                     case 'difference'
%                         [sigma_hist(:,:,k)] = updateCovariance(sigma_hist(:,:,k),dtheta(:,k),max(i,options.AM.init_memory_length),...
%                             sqrt(2)/options.AM.memory_length,options.AM.min_regularisation);
%                         mu_hist(:,k) = theta(:,k);
%                 end
%             
%                 % Proposal update
%                 if strcmp(options.proposal_scheme,'AM') && (mod(i,options.AM.adaption_interval) == 0)
%                     if acc/i < options.AM.min_acc
%                         sigma_scale(k) = sigma_scale(k)*options.AM.adap_sigma_scale;
%                     elseif acc/i > options.AM.max_acc
%                         sigma_scale(k) = sigma_scale(k)/options.AM.adap_sigma_scale;
%                     end
%                     Sigma(:,:,k) = sigma_scale(k)*sigma_hist(:,:,k);
% 
%                     % Regularisation
%                     [~,p] = cholcov(Sigma(:,:,k),0);
%                     if p ~= 0
%                         Sigma(:,:,k) = Sigma(:,:,k) + options.AM.min_regularisation*eye(parameters.number);
%                     end
%                 end
%             end
%             
%             % Store
%             if (mod(i-options.nsimu_warmup,options.thinning) == 0) && (i > options.nsimu_warmup)
%                 j = j + 1;
%                 parameters.S.PT.par(:,j,:) = theta;
%                 parameters.S.PT.logPost(j,:) = logP;
%                 
%                 str = num2str(100*acc(1)/i,'%.2f');
%                 for k = 2:options.MC.n_temps
%                     str = [str ', ' num2str(100*acc(k)/i,'%.2f')];
%                 end
%                 disp(['acc = ' str ' %']);
% 
%             end        
%         end
%         % Reduction
%         parameters.S.PT.par = parameters.S.PT.par(:,1:j,:);
%         parameters.S.PT.logPost = parameters.S.PT.logPost(1:j,:);
%         parameters.S.par = parameters.S.PT.par(:,1:j,end);
%         parameters.S.logPost = parameters.S.PT.logPost(1:j,end);
end

%% Visualization of results
if strcmp(options.mode,'visual')
    try
        % Diagnosis plots
        plotMCMCdiagnosis(parameters,'log-posterior',fh_logPost_trace);
        plotMCMCdiagnosis(parameters,'parameters',fh_par_trace);

        % Parameter distribution
        plotParameterSamples(parameters,'1D',fh_par_dis_1D,[],options.plot_options);
        %plotParameterSamples(parameters,'2D',fh_par_dis_2D,[],options.plot_options);

        % Chain statistics
        chainstats(parameters.S.par');
    end
end

%% Output
switch options.mode
    case {'visual','text'}, disp('-> Sampling FINISHED.');
    case 'silent' % no output
end

end



%% Proposal calculating function
% This function determines the mean and covariance of the MALA proposal and
% regularised / adaptive variants of it.
%   theta ... current parameter vector
%   grad ... gradient of log-posterior
%   H ... hessian or hessian approximation of log-posterior
%   beta ... regularisation parameter
%   w_hist ... weighting of history
%       = 0 ... only MALA
%       = 1 ... only adaptive Metropolis
%   mu_hist ... temporal average of theta
%   sigma_hist ... temporal covariance of theta
%   lb ... lower bound for theta
%   ub ... upper bound for theta
function [mu,Sigma] = getProposal(theta,grad,H,beta,w_hist,mu_hist,sigma_hist,lb,ub)

% Planning of MALA step
if w_hist ~= 1
    % Regularisation
    [~,p] = cholcov(-H,0);
    if p ~= 0
        k = 0;
        while p ~= 0
            H_k = H - 10^k*beta*eye(length(theta));
            [~,p] = cholcov(-H_k,0);
            k = k+1;
        end
        H = H_k;
    end
    
    % Newton step
    Sigma_MALA = -inv(H);
    Sigma_MALA = 0.5*(Sigma_MALA+Sigma_MALA');
    mu_MALA = theta - H\grad;
end

% Interpolation between
% a)  MALA (w_hist = 0) and
% b)  adaptive Metropolis with stabilized mean (w_hist = 1) 
if w_hist == 0
    % a) MALA
    Sigma = Sigma_MALA;
    mu = mu_MALA;
elseif w_hist == 1
    % b) AM
    Sigma = sigma_hist;
    mu = mu_hist;
else
    % c) Hybrid of MALA and AM
    Sigma = inv((1-w_hist)*inv(Sigma_MALA) + w_hist*inv(sigma_hist));
    Sigma = 0.5*(Sigma+Sigma');
    mu = Sigma*((1-w_hist)*inv(Sigma_MALA)*mu_MALA + w_hist*inv(sigma_hist)*mu_hist);
end

end


%% Objetive function interface
% This function is used as interface to the user-provided objective
% function. It adapts the sign and supplies the correct number of outputs.
% Furthermore, it catches errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = logPost(theta,fun,type,sign)

switch sign
    case 'negative'
        s = -1;
    case 'positive'
        s = +1;
end

try
    switch nargout
        case 1
            J = fun(theta);
            if isnan(J)
                error('J is NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J};
                case 'negative log-posterior' , varargout = {s*-J};
            end
        case 2
            [J,G] = fun(theta);
            if max(isnan([J;G(:)]))
                error('J and/or G contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:)};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:)};
            end
        case 3
            [J,G,H] = fun(theta);
            if max(isnan([J;G(:);H(:)]))
                error('J, G and/or H contain a NaN.');
            end
            switch type
                case 'log-posterior'          , varargout = {s* J,s* G(:),s* H};
                case 'negative log-posterior' , varargout = {s*-J,s*-G(:),s*-H};
            end
    end
catch error_msg
    disp(['Objective function evaluation failed because: ' error_msg.message]);
    switch nargout
        case 1
            varargout = {-s*inf};
        case 2
            varargout = {-s*inf,zeros(length(theta),1)};
        case 3
            varargout = {-s*inf,zeros(length(theta),1),zeros(length(theta))};
    end
end

end

%% Objetive function interface for parallel tempering
% This function is used as interface to the user-provided objective
% function for the parallel tempering sampling algorithm. It adapts the 
% sign and supplies the correct number of outputs. Furthermore, it catches
% errors in the user-supplied objective function.
%   theta ... parameter vector
%   fun ... user-supplied objective function
%   type ... type of user-supplied objective function
function varargout = logPost_pt(theta,fun,type,sign)

switch sign
    case 'negative'
        s = -1;
    case 'positive'
        s = +1;
end

try
    switch nargout
        case 2
            [logL,logPrior] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {s* logL,s* logPrior};
                case 'negative log-posterior' , varargout = {s*-logL,s*-logPrior};
            end
        case 6
            [logL,logPrior,dlogL,dlogPrior,ddlogL,ddlogPrior] = fun(theta);
            switch type
                case 'log-posterior'          , varargout = {s* logL,s* logPrior,s* dlogL,s* dlogPrior,s* ddlogL,s* ddlogPrior};
                case 'negative log-posterior' , varargout = {s*-logL,s*-logPrior,s*-dlogL,s*-dlogPrior,s*-ddlogL,s*-ddlogPrior};
            end
    end
catch error_msg
    disp(['Objective function evaluation failed because: ' error_msg.message]);
    switch nargout
        case 2
            varargout = {-s*inf,-s*inf};
        case 6
            varargout = {-s*inf,-s*inf,zeros(length(theta),1),zeros(length(theta),1),zeros(length(theta)),zeros(length(theta))};
    end
end

end
